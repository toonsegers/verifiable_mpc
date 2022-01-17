"""Implements https://eprint.iacr.org/2020/152.

This module implements the circuit satisfiability protocol
(Protocol 8) from ``Compressed Sigma-Protocol Theory and
Practical Application to Plug & Play Secure Algorithmics''.

Example:
    >>> proof = circuit_sat_prover(generators, code, x, ...)
    >>> tests = circuit_sat_verifier(proof, ...)

"""

import os
import sys
from enum import Enum
import logging
from random import SystemRandom
import re

import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.compressed_pivot as compressed_pivot
import verifiable_mpc.ac20.knowledge_of_exponent as koe
from verifiable_mpc.ac20.recombine import _recombination_vectors
import verifiable_mpc.tools.code_to_r1cs as c2r
import verifiable_mpc.tools.qap_creator as qc

prng = SystemRandom()

logger_cs = logging.getLogger("circuit_sat")
logger_cs.setLevel(logging.INFO)

logger_cs_hin = logging.getLogger("circuit_sat_hash_inputs")
logger_cs_hin.setLevel(logging.INFO)

logger_cs_hout = logging.getLogger("circuit_sat_hash_outputs")
logger_cs_hout.setLevel(logging.INFO)


class PivotChoice(Enum):
    """Select pivot proof system. """

    pivot = 1
    compressed = 2
    koe = 3


def create_generators(g_length, pivot_choice, group=None, progress_bar=False):
    """Create generators g, h, k.

    Ensure that, for parties that may subsequently act as provers,
    finding nontrivial linear relations between them is computationally
    as hard as computing discrete logarithms in the group.

    Keyword arguments:
    group -- the group(s) to create generators for
    progress_bar -- boolean to indicate usage of progress bar

    """
    def create_g_h():
        # Parameter `group` is required for (compressed) pivot 
        assert group is not None
        h = group.generator
        random_exponents = []
        random_exponents = list(prng.randrange(1, group.order) for i in range(g_length))

        # g = [h**i for i in random_exponents]
        g = []
        n = len(random_exponents)
        for i, r in enumerate(random_exponents):
            g.append(h ** r)
            if progress_bar:
                print(f"Generating keys: {round(100*i/(n+1))}%", end="\r")            

        return g, h

    if pivot_choice == PivotChoice.pivot:
        g, h = create_g_h()
        generators = {"g": g, "h": h}
    elif pivot_choice == PivotChoice.compressed:
        g, h = create_g_h()
        k = h ** prng.randrange(1, group.order)
        generators = {"g": g, "h": h, "k": k}
    elif pivot_choice == PivotChoice.koe and isinstance(group, list):
        group1 = group[0]
        group2 = group[1]
        order = group1.order
        _g1 = group1.generator  # BN256 regular
        _g2 = group2.generator  # BN256 twist 
        generators = koe.trusted_setup(_g1, _g2, g_length, order, progress_bar)
    else:
        raise NotImplementedError

    return generators


def input_length_power_of_2(x, code, pad_with=0):
    """Pad code and x to meet requirements of compressed pivot protocol.

    Arguments:
    x -- intput vector
    code -- string with Python function
    """

    # Calculate length of z-vector in circuit sat protocol
    inputs, body = c2r.extract_inputs_and_body(c2r.parse(code))
    flatcode = c2r.flatten_body(body)
    mul_indices_of_flatcode = mul_in_flatcode(flatcode)
    m = len(mul_indices_of_flatcode)
    z_len = len(inputs) + 3 + 2 * m

    # Add padding to x such that len(z) + 1 is a power of 2
    if not bin(z_len + 1).count("1") == 1:
        padding = next_power_of_2(z_len) - z_len - 1
    else:
        padding = 0
    padded_x = x + [type(x[0])(pad_with)] * padding
    new_z_len = len(padded_x) + 3 + 2 * m
    assert (
        bin(new_z_len + 1).count("1") == 1
    ), "This implementation requires n+1 to be power of 2 (else, use padding with zeros)."

    # Add padding to function definition in code
    def_line = [item for item in code.split("\n") if item.startswith("def")][0]
    function_params = def_line[def_line.find("(") + 1: def_line.find(")")]
    input_vars = [x.strip() for x in function_params.split(",")]
    new_input_vars = input_vars + ["padding_" + str(i) for i in range(padding)]
    padded_code = code.replace(function_params, ", ".join(new_input_vars),1)
    return padded_x, padded_code, new_z_len


def calculate_witness(code, input_vars):
    inputs, body = c2r.extract_inputs_and_body(c2r.parse(code))
    flatcode = c2r.flatten_body(body)
    witness = c2r.assign_variables(inputs, input_vars, flatcode)
    return witness


def mul_in_flatcode(flatcode):
    return [i for i, line in enumerate(flatcode) if line[0] == "*"]


def express_as_x_or_gamma(symbol, flatcode, varnames, n):
    mul_indices_of_flatcode = mul_in_flatcode(flatcode)
    m = len(mul_indices_of_flatcode)
    symbols_for_x = [sym for sym in varnames[1 : n + 1]]
    symbols_for_gammas = [flatcode[i][1] for i in mul_indices_of_flatcode]

    def get_coeff(input_string):
        coeff = re.findall("^\d*", str(input_string))[0]
        if coeff == "":
            coeff = 1
        else:
            coeff = int(coeff)
        return coeff

    def get_symbol(input_string):
        symbol = re.sub("^\d*", "", str(input_string))
        if symbol == "":
            symbol = 1
        return symbol

    def corresponding_coeff(symbol):
        return coeffs[symbols.index(symbol)]

    def replace_by_form(symbols, coeffs, varnames):
        z_symbols = symbols_for_x + ["f0", "g0", "h0"] + symbols_for_gammas
        form = pivot.LinearForm(
            [corresponding_coeff(z_i) if z_i in symbols else 0 for z_i in z_symbols]
            + [0] * m
        ) + sum(corresponding_coeff(s) for s in symbols if isinstance(s, int))
        return form

    def split(symbol):
        symbol_indices_of_flatcode = [
            i for i, line in enumerate(flatcode) if line[1] == symbol
        ]

        if len(symbol_indices_of_flatcode) > 1:
            raise NotImplementedError(
                "Code assumes that symbol only gets assigned once."
            )
        elif len(symbol_indices_of_flatcode) == 0:
            raise ValueError("Symbol does not get assigned a value in flatcode.")
        else:
            pass
        symbol_ix = symbol_indices_of_flatcode[0]

        if flatcode[symbol_ix][0] == "+":
            sym_lhs = get_symbol(flatcode[symbol_ix][2])
            sym_rhs = get_symbol(flatcode[symbol_ix][3])
            coeff_lhs = get_coeff(flatcode[symbol_ix][2])
            coeff_rhs = get_coeff(flatcode[symbol_ix][3])
        elif flatcode[symbol_ix][0] == "-":
            sym_lhs = get_symbol(flatcode[symbol_ix][2])
            sym_rhs = get_symbol(flatcode[symbol_ix][3])
            coeff_lhs = get_coeff(flatcode[symbol_ix][2])
            coeff_rhs = -1 * get_coeff(flatcode[symbol_ix][3])
        elif flatcode[symbol_ix][0] == "set":
            sym_lhs = get_symbol(flatcode[symbol_ix][2])
            sym_rhs = None
            coeff_lhs = get_coeff(flatcode[symbol_ix][2])
            coeff_rhs = None
        else:
            raise NotImplementedError

        # Code requires the z vector to be of form [x1, x2, .., f(0), g(0), h(0) = gamma1, h(1), gamma2, ...]
        symbols = []
        coeffs = []
        for symbol, coeff in [(sym_lhs, coeff_lhs), (sym_rhs, coeff_rhs)]:
            if symbol is None:
                pass
            elif symbol in symbols_for_x:
                symbols.append(symbol)
                coeffs.append(coeff)
            elif symbol in symbols_for_gammas:
                symbols.append(symbol)
                coeffs.append(coeff)
            elif isinstance(symbol, int):
                symbols.append(symbol)
                coeffs.append(coeff)
            elif not type(symbol) in [int, str]:
                raise NotImplementedError(
                    "Not able to parse other value types than ints and strings. (e.g. not able to parse sectypes)"
                )
            else:
                symbols_from_recursion, coeffs_from_recursion = split(symbol)
                symbols.extend(symbols_from_recursion)
                coeffs.extend(coeffs_from_recursion)
        return symbols, coeffs

    symbols, coeffs = split(symbol)
    form = replace_by_form(symbols, coeffs, varnames)
    return form


def mul_gates_for_splitting(flatcode, varnames, n):
    """Checks all mul-gates of flatcode for terms that require splitting;

    Returns symbol and index of terms that are not expressed as terms of z (x or gamma terms) or of type integer.
    """

    mul_indices_of_flatcode = mul_in_flatcode(flatcode)
    symbols_for_x = [sym for sym in varnames[1: n + 1]]
    symbols_for_gammas = [flatcode[i][1] for i in mul_indices_of_flatcode]
    z_symbols = symbols_for_x + symbols_for_gammas

    requires_splitting = []
    for j in [2, 3]:
        for ix in mul_indices_of_flatcode:
            symbol = flatcode[ix][j]
            if not (symbol in z_symbols or isinstance(symbol, int)):
                if j == 2:
                    wiretype = "left"
                elif j == 3:
                    wiretype = "right"
                else:
                    raise ValueError
                requires_splitting.append((symbol, ix, wiretype))

    return requires_splitting


def lagrange(gf, lagr_range, c):
    return _recombination_vectors(gf, lagr_range, (c,))[0]


def create_fgh_linear_forms(
    r1cs, c, varnames, flatcode, mul_indices_of_flatcode, n, m, gf
):
    """Create linear forms corresponding to f, g and h

    """

    def create_linear_form(M, wiretype):
        """Return linear form corresponding to poly f, g or h

        Input: Matrix of left, right or output wires, A, B or C. (represented by M)
        Output: linear form
        Code requires the varnames/M[ix] vector to be of form [1, *inputs, output, *dummy variables]
        """
        if wiretype == "left":
            poly_at_0_index = n
            lagr_range = range(m + 1)
        elif wiretype == "right":
            poly_at_0_index = n + 1
            lagr_range = range(m + 1)
        elif wiretype == "out":
            poly_at_0_index = n + 2
            lagr_range = range(2 * m + 1)
        else:
            raise ValueError(
                f"Wiretype {wiretype} should be either 'left', 'right' or 'out'."
            )

        x_terms_j = lambda ix: pivot.LinearForm(
            [gf(M[ix][i + 1]) for i in range(n)] + [gf(0)] * (3 + 2 * m)
        )
        gamma_terms_j = lambda ix: pivot.LinearForm(
            [gf(0)] * (n + 3)
            + [gf(M[ix][gamma_ix]) for gamma_ix in gamma_indices_of_varnames]
            + [0] * m
        )

        # The two lines above miss mul-gate terms that are not x or gamma; the code block below adds those terms
        symbols_to_split = mul_gates_for_splitting(flatcode, varnames, n)
        split_terms = [
            express_as_x_or_gamma(s_tuple[0], flatcode, varnames, n)
            for s_tuple in symbols_to_split
        ]
        other_terms_j = lambda ix: sum(
            gf(M[ix][varnames.index(s_tuple[0])]) * split_terms[i]
            for i, s_tuple in enumerate(symbols_to_split)
            if (s_tuple[1] == ix and s_tuple[2] == wiretype)
        )

        uvw_form_j = (
            lambda ix: gf(M[ix][0])
            + x_terms_j(ix)
            + gamma_terms_j(ix)
            + other_terms_j(ix)
        )
        poly_at_0 = [0] * (n + 3 + 2 * m)
        poly_at_0[poly_at_0_index] = 1
        linform_0_to_m = pivot.LinearForm(poly_at_0) * lagrange(gf, lagr_range, c)[0] + sum(
            uvw_form_j(ix) * lagrange(gf, lagr_range, c)[j + 1] for j, ix in enumerate(mul_indices_of_flatcode)
        )

        if wiretype in ["left", "right"]:
            linform = linform_0_to_m
        if wiretype == "out":
            linform_0_to_2m = (
                linform_0_to_m
                + pivot.LinearForm(
                    [0] * (n + 3 + m)
                    + [1 * lagrange(gf, lagr_range, c)[1 + m + i] for i in range(m)]
                )
                + linform_0_to_m.constant
            )
            linform = linform_0_to_2m

        return linform

    A, B, C = r1cs
    gamma_indices_of_varnames = [
        varnames.index(flatcode[i][1]) for i in mul_indices_of_flatcode
    ]
    linform_f = create_linear_form(A, "left")
    linform_g = create_linear_form(B, "right")
    linform_h = create_linear_form(C, "out")

    return linform_f, linform_g, linform_h


def code_to_flatcode_and_r1cs(code):
    inputs, body = c2r.extract_inputs_and_body(c2r.parse(code))
    flatcode = c2r.flatten_body(body)
    varnames = c2r.get_var_placement(inputs, flatcode)
    r1cs = c2r.flatcode_to_r1cs(inputs, flatcode)
    return flatcode, inputs, varnames, r1cs


def _inner_prod_asymmetric(v1, v2):
    out = type(v2[-1])(0)
    for k in [i for i, e in enumerate(v1) if e != 0]:
        out += v1[k] * v2[k]
    return out


def calculate_ab_vectors(r1cs, xc, mul_indices_of_flatcode):
    A, B, C = r1cs

    A_mul_rows = [A[j] for j in mul_indices_of_flatcode]
    B_mul_rows = [B[j] for j in mul_indices_of_flatcode]
    a = [_inner_prod_asymmetric(A_j, xc) for A_j in A_mul_rows]
    b = [_inner_prod_asymmetric(B_j, xc) for B_j in B_mul_rows]

    return a, b


def calculate_fgh_polys(a, b, c, gf):  # TODO: don't pass c, not needed
    # Calculate random polynomials f, g, h
    r_a = prng.randrange(1, gf.order)
    r_b = prng.randrange(1, gf.order)
    f_poly = qc.Poly(qc.lagrange_interp_ff(a + [r_a], gf))  # TODO: check: why not random element at position 0? why at end? like on p.21 on AC20 paper, it says f(i)=alpha_i for i=1,..,m
    g_poly = qc.Poly(qc.lagrange_interp_ff(b + [r_b], gf))
    h_poly = f_poly * g_poly
    # assert c == [h_poly.eval(i+1) for i in range(m)], "Evaluations of h at 1..m not equal to vector c"
    return f_poly, g_poly, h_poly


def next_power_of_2(x):
    return 1 << (x).bit_length()


def protocol_8_excl_pivot_prover(generators, code, x, gf, use_koe=False):
    """Non-interactive implementation of Protocol 8,
    including nullity protocol, but excluding the call
    to the compressed pivot protocol.

    """
    if "g" in generators:
        g = generators["g"]
        h = generators["h"]
    elif "pp_lhs" in generators:
        use_koe = True
        pp = generators
    else:
        raise NotImplementedError

    n = len(x)
    xc = calculate_witness(code, x)
    proof = {}
    flatcode, inputs, varnames, r1cs = code_to_flatcode_and_r1cs(code)

    mul_indices_of_flatcode = mul_in_flatcode(flatcode)
    m = len(mul_indices_of_flatcode)
    output_variables = [symbol for symbol in varnames if symbol.startswith("~out")]

    # Calculate a, b, c vectors, only for mul gates
    a, b = calculate_ab_vectors(r1cs, xc, mul_indices_of_flatcode)
    c = [a_i * b_i for a_i, b_i in zip(a, b)]

    # Calculate random polynomials f, g, h
    f_poly, g_poly, h_poly = calculate_fgh_polys(a, b, c, gf)

    # Prover commits to (x, f(0), g(0), h(0), h(1), .., h(2m))
    h_evaluations = [h_poly.eval(i + 1) for i in range(2 * m)]
    z = x + [f_poly.eval(0), g_poly.eval(0), h_poly.eval(0)] + h_evaluations

    gamma = prng.randrange(1, gf.order)

    if use_koe:
        S = range(len(z))
        z_commitment_P, z_commitment_pi = koe.restriction_argument_prover(
            S, z, gamma, pp
        )
        z_commitment = {"P": z_commitment_P, "pi": z_commitment_pi}
        proof["z_commitment"] = z_commitment
    else:
        logger_cs.debug("Calculate [Z].")
        z_commitment = pivot.vector_commitment(z, gamma, g, h)
        proof["z_commitment"] = z_commitment

    logger_cs.debug("Apply first hash of circuit satisfiability protocol.")
    input_list = [z_commitment, code, "First hash circuit satisfiability protocol"]
    logger_cs_hin.debug(f"Method protocol_8_excl_pivot_prover: input_list=\n{input_list}")
    c = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_hout.debug(f"After hash, hash=\n{c}")

    # Prover's response
    y1 = f_poly.eval(c)
    y2 = g_poly.eval(c)
    y3 = h_poly.eval(c)
    assert y3 == y1 * y2

    """Create linear forms corresponding to f(c), g(c), h(c) for Pi_Nullity step
    in Protocol 6 (Protocol 8 in updated version)
    """
    linform_f, linform_g, linform_h = create_fgh_linear_forms(
        r1cs, c, varnames, flatcode, mul_indices_of_flatcode, n, m, gf
    )

    y1 = linform_f(z)
    y2 = linform_g(z)
    y3 = linform_h(z)
    assert y1 * y2 == y3
    proof["y1"] = y1
    proof["y2"] = y2
    proof["y3"] = y3

    circuits = []
    outputs = []
    for output_var in output_variables:
        circuit = express_as_x_or_gamma(output_var, flatcode, varnames, n)
        y = circuit(z)
        assert (
            y == xc[varnames.index(output_var)]
        ), f"Output of circuit {y} not equal to ~out in witness."
        circuits.append(circuit)
        outputs.append(y)
    proof["outputs"] = outputs

    lin_forms = [circuit - y for circuit, y in zip(circuits, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    # Follow prove_nullity_compressed but add more inputs to Fiat-Shamir hash
    logger_cs.debug("Apply second hash of circuit satisfiability protocol.")
    input_list = [
        y1,
        y2,
        y3,
        z_commitment,
        outputs,
        circuits,
        lin_forms,
        "Second hash circuit satisfiability protocol",
    ]
    logger_cs_hin.debug(f"Method: protocol_8_excl_pivot_prover, input_list=\n{input_list}")
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    proof["L"] = L
    return proof, z_commitment, L, z, gamma


def protocol_8_excl_pivot_verifier(proof, code, gf, use_koe=False):
    verification = {}
    y1 = proof["y1"]
    y2 = proof["y2"]
    y3 = proof["y3"]
    if not y1 * y2 == y3:
        verification["y1*y2=y3"] = False
        return verification
    else:
        verification["y1*y2=y3"] = True

    flatcode, inputs, varnames, r1cs = code_to_flatcode_and_r1cs(code)
    n = len(inputs)

    mul_indices_of_flatcode = mul_in_flatcode(flatcode)
    m = len(mul_indices_of_flatcode)
    output_variables = [symbol for symbol in varnames if symbol.startswith("~out")]

    if "P" in proof:
        use_koe = True
    if use_koe:
        z_commitment_P = proof["z_commitment"]["P"]
        z_commitment_pi = proof["z_commitment"]["pi"]
        logger_cs.debug("Apply first hash of circuit satisfiability protocol.")
        input_list = [
            z_commitment_P,
            z_commitment_pi,
            code,
            "First hash circuit satisfiability protocol",
        ]
        logger_cs_hin.debug(f"Method: protocol_8_excl_pivot_verifier (1), input_list=\n{input_list}")
        c = pivot.fiat_shamir_hash(input_list, gf.order)
        logger_cs_hout.debug(f"After hash, hash=\n{c}")
    else:
        z_commitment = proof["z_commitment"]
        logger_cs.debug("Apply first hash of circuit satisfiability protocol.")
        input_list = [z_commitment, code, "First hash circuit satisfiability protocol"]
        logger_cs_hin.debug(f"Method: protocol_8_excl_pivot_verifier (1), input_list=\n{input_list}")
        c = pivot.fiat_shamir_hash(input_list, gf.order)
        logger_cs_hout.debug(f"After hash, hash=\n{c}")

    linform_f, linform_g, linform_h = create_fgh_linear_forms(
        r1cs, c, varnames, flatcode, mul_indices_of_flatcode, n, m, gf
    )

    outputs = proof["outputs"]
    circuits = []
    for output_var in output_variables:
        circuit = express_as_x_or_gamma(output_var, flatcode, varnames, n)
        circuits.append(circuit)

    lin_forms = [circuit - output for circuit, output in zip(circuits, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    logger_cs.debug("Apply second hash of circuit satisfiability protocol.")
    input_list = [
        y1,
        y2,
        y3,
        z_commitment,
        outputs,
        circuits,
        lin_forms,
        "Second hash circuit satisfiability protocol",
    ]
    logger_cs_hin.debug(f"Method: protocol_8_excl_pivot_verifier (2), input_list=\n{input_list}")
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))

    if not L == proof["L"]:
        verification["L_wellformed_from_Cfgh_forms"] = False
        return verification, L
    else:
        verification["L_wellformed_from_Cfgh_forms"] = True

    return verification, L


def circuit_sat_prover(generators, code, x, gf, pivot_choice=PivotChoice.compressed):
    """Non-interactive implementation of Protocol 8, prover-side,
    including Nullity using arbitrary pivot.
    """
    proof, z_commitment, L, z, gamma = protocol_8_excl_pivot_prover(
        generators, code, x, gf
    )

    if pivot_choice == PivotChoice.compressed:
        pivot_proof = compressed_pivot.protocol_5_prover(
            generators, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == PivotChoice.pivot:
        g = generators["g"]
        h = generators["h"]
        pivot_proof = pivot.prove_linear_form_eval(
            g, h, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == PivotChoice.koe:
        L = proof["L"]
        P = proof["z_commitment"]["P"]
        pi = proof["z_commitment"]["pi"]
        pivot_proof, u = koe.opening_linear_form_prover(L, z, gamma, generators, P, pi)
    else:
        raise NotImplementedError
    proof["pivot_proof"] = pivot_proof

    return proof


def circuit_sat_verifier(
    proof, generators, code, gf, pivot_choice=PivotChoice.compressed
):
    """Non-interactive implementation of Protocol 8, verifier-side,
    including Nullity using arbitrary pivot.
    """
    verification, L = protocol_8_excl_pivot_verifier(proof, code, gf)

    # Run pivot protocol based on enum `pivot_choice`
    if pivot_choice == PivotChoice.compressed:
        z_commitment = proof["z_commitment"]
        pivot_proof = proof["pivot_proof"]
        pivot_verification = compressed_pivot.protocol_5_verifier(
            generators, z_commitment, L, 0, pivot_proof, gf
        )
    elif pivot_choice == PivotChoice.pivot:
        z_commitment = proof["z_commitment"]
        g = generators["g"]
        h = generators["h"]
        z, phi, c = proof["pivot_proof"]
        pivot_verification = pivot.verify_linear_form_proof(
            g, h, z_commitment, L, 0, z, phi, c
        )
    elif pivot_choice == PivotChoice.koe:
        koe_pivot_proof = proof["pivot_proof"]
        pivot_verification = koe.opening_linear_form_verifier(
            L, generators, koe_pivot_proof, 0
        )
        verification["pivot_verification"] = pivot_verification
    else:
        raise NotImplementedError
    verification["pivot_verification"] = pivot_verification

    return verification


# TODO-list
# TODO-1: Avoid double work: compute z_commitment once, don't compute P (equal to z_commitment) again
# TODO-1: Apply DRY to circuit_sat sub-methods; remove redundant code
# TODO-1: Rename `symbols` to `vars`
# TODO-1: Clean up create_fgh_linear_forms
# TODO-1: Add "tag" to hashes
# TODO-3: Reduce mul-gates: represent multiply-by-constant lines in flatcode by constant in R1CS matrices
