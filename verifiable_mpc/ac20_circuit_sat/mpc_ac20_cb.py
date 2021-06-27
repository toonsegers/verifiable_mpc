""" MPC interface for circuit_sat and compressed_pivot modules.

Circuit_sat implementation of https://eprint.iacr.org/2020/152
``Compressed Sigma-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

"""

import pprint

import verifiable_mpc.ac20_circuit_sat.pivot as pivot
import verifiable_mpc.ac20_circuit_sat.circuit_builder as cb
import verifiable_mpc.ac20_circuit_sat.circuit_sat_r1cs as cs
import verifiable_mpc.ac20_circuit_sat.knowledge_of_exponent as koe
from verifiable_mpc.ac20_circuit_sat.pivot import _int as _int
from mpyc.runtime import mpc, logging
from mpyc.random import randrange
from sec_groups.fingroups import EllipticCurveElement
from sec_groups.secgroups import secure_repeat_public_base_public_output as secure_repeat
import verifiable_mpc.tools.qap_creator as qc


pp = pprint.PrettyPrinter(indent=4)

logger_cs_mpc_cb = logging.getLogger("cs_mpc")
logger_cs_mpc_cb.setLevel(logging.INFO)

logger_cs_mpc_cb_hin = logging.getLogger("cs_mpc_hash_inputs")
logger_cs_mpc_cb_hin.setLevel(logging.INFO)

logger_cs_mpc_cb_hout = logging.getLogger("cs_mpc_hash_outputs")
logger_cs_mpc_cb_hout.setLevel(logging.INFO)


def vector_commitment(x, gamma, g, h):
    """Pedersen vector commitment. Definition 1 of AC20.

    Uses secure groups to locally compute Prod_j g_j^{share_i of x_j},
    reduce result and publish this to other parties for recombination.
    """
    c = secure_repeat(g + [h], x + [gamma])
    return c


async def create_generators(group, sectype, input_length):
    # TODO: consider setup to exclude trapdoors
    h = group.generator
    random_exponents = [mpc._random(sectype) for i in range(input_length + 1)]
    kg = await mpc.gather([secure_repeat(h, u) for u in random_exponents])
    generators = {"g": kg[1:], "h": h, "k": kg[0]}
    return generators


async def protocol_4_prover(g_hat, k, Q, L_tilde, z_hat, gf, proof={}, round_i=0):
    """Non-interactive version of Protocol 4 from Section 4.2: Prover's side"""

    # Step 5: Prover calculates A, B
    half = len(g_hat) // 2
    g_hat_l = g_hat[:half]
    g_hat_r = g_hat[half:]
    z_hat_l = z_hat[:half]
    z_hat_r = z_hat[half:]
    logger_cs_mpc_cb.debug("Calculate A_i, B_i.")
    A = await vector_commitment(
        z_hat_l, _int(L_tilde([0] * half + z_hat_l)), g_hat_r, k
    )
    B = await vector_commitment(
        z_hat_r, _int(L_tilde(z_hat_r + [0] * half)), g_hat_l, k
    )
    logger_cs_mpc_cb.debug(f"Provers opened A{str(round_i)}, B{str(round_i)}")

    proof["A" + str(round_i)] = A
    proof["B" + str(round_i)] = B

    # Step 6: Prover calculates challenge
    order = k.order
    # Hash A, B and all public elements of Protocol 4

    # input_list = [A, B, g_hat, k, Q, L_tilde]
    if isinstance(A, EllipticCurveElement):
        input_list = [A.to_affine(), B.to_affine(), g_hat, k, Q.to_affine(), L_tilde]
    else:
        input_list = [A, B, g_hat, k, Q, L_tilde]

    logger_cs_mpc_cb_hin.debug(
        f"Method protocol_4_prover: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c = pivot.fiat_shamir_hash(input_list, order)
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{c}")

    # Step 7: Prover calculates following public values
    logger_cs_mpc_cb.debug("Calculate g_prime.")
    g_prime = [(g_hat_l[i] ** c) * g_hat_r[i] for i in range(half)]
    logger_cs_mpc_cb.debug("Calculate Q_prime.")
    Q_prime = A * (Q ** c) * (B ** (c ** 2))

    assert (
        L_tilde.constant == 0
    ), "Next line assumes L_tilde is a linear form, not affine form."
    c_L_tilde_l_coeffs = [coeff * gf(c) for coeff in L_tilde.coeffs[:half]]
    L_prime = pivot.LinearForm(c_L_tilde_l_coeffs) + pivot.LinearForm(
        L_tilde.coeffs[half:]
    )

    # Step 8: Prover calculates z_prime and tests if recursion is required (if z' in Z or Z^2)
    z_prime = [z_hat_l[i] + c * z_hat_r[i] for i in range(half)]
    if len(z_prime) <= 2:
        z_prime = await mpc.output(
            z_prime
        )  # TODO: perhaps not necessary to open z_prime; double check
        logger_cs_mpc_cb.debug(f"Provers opened z_prime")
        proof["z_prime"] = z_prime
        return proof
    else:
        # Step 9: Prover sends z_prime and starts recursion of protocol 4
        round_i += 1
        proof = await protocol_4_prover(
            g_prime, k, Q_prime, L_prime, z_prime, gf, proof, round_i
        )
        return proof


async def protocol_5_prover(generators, P, L, y, x, gamma, gf):
    secfld = type(x[0])
    g = generators["g"]
    h = generators["h"]
    k = generators["k"]

    proof = {}

    # Public parameters
    n = len(x)
    L, y = pivot.affine_to_linear(L, y, n)
    L.constant = await mpc.output(L.constant)

    # Open P, y, L to ensure consistent output of Fiat-Shamir hash
    y = await mpc.output(y)
    assert (
        bin(n + 1).count("1") == 1
    ), "This implementation requires n+1 to be power of 2 (else, use padding with zeros)."
    logger_cs_mpc_cb.debug(f"Provers opened y.")

    # Non-interactive proof
    # Step 1: Prover calculates ("announcement") t, A
    order = gf.order
    r = list(mpc._random(secfld) for i in range(n))
    rho = mpc._random(secfld)
    t = L(r)
    logger_cs_mpc_cb.debug("Calculate A.")
    A = await vector_commitment(r, rho, g, h)
    t = await mpc.output(t)
    logger_cs_mpc_cb.debug(f"Provers opened t, A.")
    proof["t"] = t
    proof["A"] = A

    # Step 2: Prover computes challenge
    # input_list = [t, A, generators, P, L, y]
    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.to_affine(), generators, P.to_affine(), L, y]
    else:
        input_list = [t, A, generators, P, L, y]

    logger_cs_mpc_cb_hin.debug(
        f"Method protocol_5_prover: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c0 = pivot.fiat_shamir_hash(
        input_list + [0] + ["First hash of compressed pivot"], order
    )
    c1 = pivot.fiat_shamir_hash(
        input_list + [1] + ["First hash of compressed pivot"], order
    )
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{c0}, {c1}")

    # Step 3: Prover calculates
    z = [c0 * x_i + r[i] for i, x_i in enumerate(x)]
    phi = c0 * gamma + rho
    z_hat = z + [phi]

    # Step 4: Prover calculates following public variables
    g_hat = g + [h]
    logger_cs_mpc_cb.debug("Calculate Q.")
    Q = A * (P ** c0) * (k ** _int(c1 * (c0 * y + t)))
    L_tilde = pivot.LinearForm(L.coeffs + [0]) * c1
    # assert L(z)*c1 == L_tilde(z_hat)
    proof = await protocol_4_prover(g_hat, k, Q, L_tilde, z_hat, gf, proof)
    return proof


def calculate_fgh_polys(a, b, c, gf, secfld):
    # Calculate random polynomials f, g, h
    r_a = mpc._random(secfld)
    r_b = mpc._random(secfld)
    logger_cs_mpc_cb.debug("Calculate f_poly.")
    f_poly = qc.Poly(qc.lagrange_interp_ff(a + [r_a], gf))
    logger_cs_mpc_cb.debug("Calculate g_poly.")
    g_poly = qc.Poly(qc.lagrange_interp_ff(b + [r_b], gf))
    logger_cs_mpc_cb.debug("Calculate h_poly.")
    h_poly = f_poly * g_poly
    logger_cs_mpc_cb.debug("Done calculating f, g, h.")
    # assert c == [h_poly.eval(i+1) for i in range(m)], "Evaluations of h at 1..m not equal to vector c"
    return f_poly, g_poly, h_poly


# START########################## from mpyc.thresha -> modified _recombination_vector ##########
import functools


@functools.lru_cache(maxsize=None)
def _recombination_vectors(field, xs, xr):
    """Compute and store recombination vectors.

    Recombination vectors depend on the field, the x-coordinates xs
    of the shares and the x-coordinates xr of the recombination points.
    """
    modulus = field.modulus
    xs = [x % modulus for x in xs]  # also for conversion from
    xr = [x % modulus for x in xr]  # int to type(modulus)
    d = [None] * len(xs)
    for i, x_i in enumerate(xs):
        q = 1
        for j, x_j in enumerate(xs):
            if i != j:
                q *= x_i - x_j
                q %= modulus
        d[i] = q
    matrix = [None] * len(xr)
    for r, x_r in enumerate(xr):
        matrix[r] = [None] * len(xs)
        p = 1
        for j, x_j in enumerate(xs):
            p *= x_r - x_j
            p %= modulus
        p = field(p)
        for i, x_i in enumerate(xs):
            matrix[r][i] = (p / field((x_r - x_i) * d[i])).value
    return matrix


def recombine(field, points, x_rs=0):  ##### ONLY for shares that are single numbers
    """Recombine shares given by points into secrets.

    Recombination is done for x-coordinates x_rs.
    """
    xs, shares = list(zip(*points))
    if not isinstance(x_rs, list):
        x_rs = (x_rs,)
    m = len(shares)
    width = len(x_rs)
    T_is_field = isinstance(shares[0], field)  # all elts assumed of same type
    vector = _recombination_vectors(field, xs, tuple(x_rs))
    sums = [0] * width
    for i in range(m):
        s = shares[i]
        if T_is_field:
            s = s.value
        # type(s) is int or gfpx.Polynomial
        for r in range(width):
            sums[r] += s * vector[r][i]
    for r in range(width):
        sums[r] = field(sums[r])
    if isinstance(x_rs, tuple):
        return sums[0]

    return sums


# END########################## from mpyc.thresha -> modified _recombination_vector ##########


async def protocol_8_excl_pivot_prover(generators, circuit, x, gf, use_koe=False):
    """Non-interactive implementation of Protocol 8,
    including nullity protocol, but excluding the call
    to the compressed pivot protocol.

    """
    secfld = type(x[0])
    if "g" in generators:
        g = generators["g"]
        h = generators["h"]
    elif "pp_lhs" in generators:
        use_koe = True
        pp = generators
    else:
        raise NotImplementedError

    logger_cs_mpc_cb.debug(f"Calculate witness.")
    n = len(x)
    assert n == circuit.input_ct
    proof = {}
    m = circuit.mul_ct

    # Calculate a, b, c vectors, only for mul gates
    logger_cs_mpc_cb.debug(f"Calculate a, b, c vectors.")
    a, b, c = circuit.multiplication_triples(x)

    logger_cs_mpc_cb.debug(f"Calculate z.")
    f0 = mpc._random(secfld)
    g0 = mpc._random(secfld)
    a = [f0] + a
    b = [g0] + b
    logger_cs_mpc_cb.debug(f"Calculate z: a, b = await mpc.gather(a,b)")
    a, b = await mpc.gather(a, b)
    logger_cs_mpc_cb.debug(f"Calculate z: fs = recombine(...)")
    fs = recombine(gf, list(zip(range(m + 1), a)), list(range(m + 1, 2 * m + 1)))
    logger_cs_mpc_cb.debug(f"Calculate z: gs = recombine(...)")
    gs = recombine(gf, list(zip(range(m + 1), b)), list(range(m + 1, 2 * m + 1)))
    logger_cs_mpc_cb.debug(f"Calculate z: hs = list(...)")
    hs = list(map(secfld, await mpc.schur_prod(fs, gs)))
    z = x + [f0, g0, f0 * g0] + c + hs

    gamma = mpc._random(secfld)

    if use_koe:
        S = range(len(z))
        z_commitment_P, z_commitment_pi = koe.restriction_argument_prover(
            S, z, gamma, pp
        )
        z_commitment = {"P": z_commitment_P, "pi": z_commitment_pi}
        proof["z_commitment"] = z_commitment
    else:
        logger_cs_mpc_cb.debug("Calculate commitment for z, denoted by [z].")
        z_commitment = await vector_commitment(z, gamma, g, h)
        # logger_cs_mpc_cb.debug("Calculated [z] in secret shared form.")
        logger_cs_mpc_cb.debug(f"Opened [z] (z_commitment).")
        proof["z_commitment"] = z_commitment

    logger_cs_mpc_cb.debug("Apply first hash of circuit satisfiability protocol.")
    input_list = [z_commitment, str(circuit), "First hash circuit satisfiability protocol"]
    logger_cs_mpc_cb_hin.debug(
        f"Method protocol_8_excl_pivot_prover (1): input_list={input_list}"
    )
    c = pivot.fiat_shamir_hash(
        input_list, gf.order
    )  # TODO: name clash c = a * b from above
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{c}")

    # Prover's response
    linform_f = cb.calculate_fg_form(circuit, wire=0, challenge=c, gf=gf)
    linform_g = cb.calculate_fg_form(circuit, wire=1, challenge=c, gf=gf)
    linform_h = cb.calculate_h_form(circuit, c, gf)

    y1 = linform_f(z)
    y2 = linform_g(z)
    y3 = linform_h(z)
    y1 = await mpc.output(y1, raw = True)  # raw = True ensures that we return field elements
    y2 = await mpc.output(y2, raw = True)
    y3 = await mpc.output(y3, raw = True)

    assert y1 * y2 == y3
    proof["y1"] = y1
    proof["y2"] = y2
    proof["y3"] = y3

    circuit_forms = cb.calculate_circuit_forms(circuit)
    circuit_forms = [cb.convert_to_ac20(f, circuit) for f in circuit_forms]
    outputs = circuit(x)
    outputs = await mpc.output(outputs)
    proof["outputs"] = outputs

    lin_forms = [form - y for form, y in zip(circuit_forms, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    # Follow prove_nullity_compressed but add more inputs to Fiat-Shamir hash
    logger_cs_mpc_cb.debug("Apply second hash of circuit satisfiability protocol.")
    input_list = [
        y1,
        y2,
        y3,
        z_commitment,
        outputs,
        circuit_forms,
        lin_forms,
        "Second hash circuit satisfiability protocol",
    ]
    logger_cs_mpc_cb_hin.debug(
        f"Method protocol_8_excl_pivot_prover (2): input_list={input_list}"
    )
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    proof["L"] = L
    return proof, z_commitment, L, z, gamma


async def prove_linear_form_eval(g, h, P, L, y, x, gamma, gf):
    """Sigma protocol Pi_s (protocol 2) from AC20.

    Non-interactive version.
    """
    secfld = type(x[0])
    n = len(x)
    L, y = pivot.affine_to_linear(L, y, n)
    y = await mpc.output(y)

    r = list(mpc._random(secfld) for i in range(n))
    rho = mpc._random(secfld)

    t = L(r)
    A = await vector_commitment(r, rho, g, h)
    t = await mpc.output(t)
    logger_cs_mpc_cb.debug(f"Provers opened t.")
    logger_cs_mpc_cb.debug(f"Provers opened A.")

    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.to_affine(), g, h, P.to_affine(), L, y]
    else:
        input_list = [t, A, g, h, P, L, y]
    logger_cs_mpc_cb_hin.debug(f"Method prove_linear_form_eval: input_list={input_list}.")
    c = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{c}")
    z = [c * x_i + r[i] for i, x_i in enumerate(x)]
    # phi = (c*gamma + rho) % gf.order
    phi = c * gamma + rho

    z = await mpc.output(z)
    phi = await mpc.output(phi)

    return (
        z,
        phi,
        c,
    )  # TODO: check if it's correct to return c as well (Required to reconstruct A in non-interactive proof)


async def circuit_sat_prover(
    generators, code, x, gf, pivot_choice=cs.PivotChoice.compressed
):
    """Non-interactive implementation of Protocol 8, prover-side,
    including Nullity using compressed pivot (Protocol 5).
    """
    logger_cs_mpc_cb.debug(f"Enter circuit_sat_prover. pivot_choice={pivot_choice}")

    logger_cs_mpc_cb.debug(f"Start protocol 8, excluding pivot proof.")
    proof, z_commitment, L, z, gamma = await protocol_8_excl_pivot_prover(
        generators, code, x, gf
    )

    if pivot_choice == cs.PivotChoice.compressed:
        pivot_proof = await protocol_5_prover(
            generators, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == cs.PivotChoice.pivot:
        g = generators["g"]
        h = generators["h"]
        pivot_proof = await prove_linear_form_eval(
            g, h, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == cs.PivotChoice.koe:
        L = proof["L"]
        P = proof["z_commitment"]["P"]
        pi = proof["z_commitment"]["pi"]
        pivot_proof, u = koe.opening_linear_form_prover(L, z, gamma, generators, P, pi)
    else:
        raise NotImplementedError
    proof["pivot_proof"] = pivot_proof

    return proof


# TODO-list
# TS-1: Include option that provers don't learn output y (or multiple outputs) and only keep y in secret shared form while proving.
# TS-1: Do y and outputs have to be part of the hash? its in second hash of protocol 8 (circuit sat) and hash of protocol 5 (compressed pivot)
# TS-1: Can the prove be forged when y is not part of the hash?
# TS-1: Avoid double work: compute z_commitment once, don't compute P (equal to z_commitment) again
# TS-1: Enable sectype operations for SecureEdwardsGroup as it does not inherit from Share but from Point
# * Also in _int() helper method?
# TS-2: In method output(), also catch SecEdwardsGroup types (check spelling)
