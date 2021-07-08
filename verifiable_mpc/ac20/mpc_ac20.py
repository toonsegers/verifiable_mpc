""" MPC interface for circuit_sat and compressed_pivot modules.

Circuit_sat implementation of https://eprint.iacr.org/2020/152
``Compressed Sigma-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

"""

import mpyc.mpctools as mpctools
from mpyc.runtime import mpc, logging
from sec_groups.fingroups import EllipticCurveElement
from sec_groups.secgroups import secure_repeat_public_base_public_output as secure_repeat
import verifiable_mpc.tools.qap_creator as qc
import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.circuit_sat_r1cs as cs
import verifiable_mpc.ac20.knowledge_of_exponent as koe
from verifiable_mpc.ac20.recombine import recombine, _recombination_vectors
from verifiable_mpc.ac20.pivot import _int as _int


logger_cs_mpc = logging.getLogger("cs_mpc")
logger_cs_mpc.setLevel(logging.INFO)

logger_cs_mpc_hin = logging.getLogger("cs_mpc_hash_inputs")
logger_cs_mpc_hin.setLevel(logging.INFO)

logger_cs_mpc_hout = logging.getLogger("cs_mpc_hash_outputs")
logger_cs_mpc_hout.setLevel(logging.INFO)


def list_mul(x):
    return mpctools.reduce(type(x[0]).operation, x)


def vector_commitment(x, gamma, g, h):
    """ Pedersen vector commitment. Definition 1 of AC20.

    Uses secure groups to locally compute Prod_j g_j^{share_i of x_j},
    reduce result and publish this to other parties for recombination.
    """
    c = secure_repeat(g + [h], x + [gamma])
    return c


async def create_generators(group, sectype, input_length):
    logger_cs_mpc.warning("Circuit reconstruction using R1CS is deprecated. Consider using circuit_builder.")
    # TODO: consider setup to exclude trapdoors
    h = group.generator
    random_exponents = [mpc._random(sectype) for i in range(input_length+1)]
    kg = await mpc.gather([secure_repeat(h, u) for u in random_exponents])
    generators = {"g": kg[1:], "h": h, "k": kg[0]}
    return generators


async def koe_trusted_setup(group, sectype, input_length, progress_bar=False):
    group1 = group[0]
    group2 = group[1]
    order = group1.order
    _g1 = group1.generator  # BN256 regular
    _g2 = group2.generator  # BN256 twist 

    g_exp = mpc._random(sectype)  # TODO: sample > 0
    alpha = mpc._random(sectype)  
    z = mpc._random(sectype)
    g1 = await secure_repeat(_g1, g_exp)
    g2 = await secure_repeat(_g2, g_exp * alpha)

    pp_lhs = []
    pp_rhs = []
    g1_base = g1
    g2_base = g2
    for i in range(2 * input_length):
        g1 = await secure_repeat(g1, z)
        g2 = await secure_repeat(g2, z)

        pp_lhs.append(g1_base)
        pp_rhs.append(g2_base)

        if progress_bar:
            print(f"Generating keys: {round(100*i/(2*input_length+1))}%", end="\r")

    pp = {"pp_lhs": pp_lhs, "pp_rhs": pp_rhs}
    return pp


async def koe_restriction_argument_prover(S, x, gamma, pp):
    """Resriction argument from [Gro10]: Prover's side.

    Extension of the Knowledge Commitment Scheme, but proofs that a subset
    of indices S of {0, ..., n-1}, corresponding to vector x, was used.
    """

    # Prover commits only to S-indices of x with (new) commitment P, with secret random gamma
    # P = (pp["pp_lhs"][0] ** gamma) * list_mul([pp["pp_lhs"][i + 1] ** x[i] for i in S])
    P = (await secure_repeat(pp["pp_lhs"][0], gamma)) * list_mul([await secure_repeat(pp["pp_lhs"][i + 1], x[i]) for i in S])

    # Trusted party or verifier samples beta; only required in designated verifier scenario
    # beta = prng.randrange(1, order)
    # sigma = (g1**beta, [(g2**(z**i))**beta for i in S])

    # Prover computes pi^alpha (slight deviation from Thomas' notes, who computes pi)
    # pi = (pp["pp_rhs"][0] ** gamma) * list_mul([pp["pp_rhs"][i + 1] ** x[i] for i in S])
    pi = (await secure_repeat(pp["pp_rhs"][0], gamma)) * list_mul([await secure_repeat(pp["pp_rhs"][i + 1], x[i]) for i in S])
    return P, pi


async def koe_opening_linear_form_prover(L, x, gamma, pp, P=None, pi=None):
    """Opening linear forms using an adaptation of the multiplication argument
    of [Gro10].

    """
    proof = {}
    n = len(x)
    S = range(n)
    # """ Run Restriction argument on (P, {1, ..., n}; x, gamma) to show that
    # (P, P_bar) is indeed a commitment to vector x in Z_q
    # """
    if P is None:
        P, pi = await koe_restriction_argument_prover(S, x, gamma, pp)
    proof["P"] = P
    proof["pi"] = pi

    # Prover computes coefficients of c_poly and Q, sends Q to verifier
    u = L(x)
    L_linear, u_linear = pivot.affine_to_linear(L, u, n)

    c_poly_lhs = qc.Poly([gamma] + [x_i for x_i in x])
    c_poly_rhs = qc.Poly([L_linear.coeffs[n - (j + 1)] for j in range(n)])
    c_poly = c_poly_lhs * c_poly_rhs

    # assert u_linear == c_poly.coeffs[n], "L(x) not equal to n-th coefficient of c_poly"
    c_bar = c_poly.coeffs
    c_bar[n] = 0  # Ensure c_bar[n] is also a sectype (for secure_repeat)
    assert len(pp["pp_lhs"]) == 2 * n
    # Q = list_mul([g_i ** (-1 * c_bar[i]) for i, g_i in enumerate(pp["pp_lhs"])])
    Q = list_mul([await secure_repeat(g_i, (c_bar[i] * -1)) for i, g_i in enumerate(pp["pp_lhs"])])
    proof["Q"] = Q
    return proof, u


async def protocol_4_prover(g_hat, k, Q, L_tilde, z_hat, gf, proof={}, round_i=0):
    """ Non-interactive version of Protocol 4 from Section 4.2: Prover's side."""

    # Step 5: Prover calculates A, B
    half = len(g_hat) // 2
    g_hat_l = g_hat[:half]
    g_hat_r = g_hat[half:]
    z_hat_l = z_hat[:half]
    z_hat_r = z_hat[half:]
    logger_cs_mpc.debug("Calculate A_i, B_i.")
    A = await vector_commitment(
        z_hat_l, _int(L_tilde([0] * half + z_hat_l)), g_hat_r, k
    )
    B = await vector_commitment(
        z_hat_r, _int(L_tilde(z_hat_r + [0] * half)), g_hat_l, k
    )
    logger_cs_mpc.debug(f"Provers opened A{str(round_i)}, B{str(round_i)}")

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

    logger_cs_mpc_hin.debug(
        f"Method protocol_4_prover: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c = pivot.fiat_shamir_hash(input_list, order)
    logger_cs_mpc_hout.debug(f"After hash, hash=\n{c}")

    # Step 7: Prover calculates following public values
    logger_cs_mpc.debug("Calculate g_prime.")
    g_prime = [(g_hat_l[i] ** c) * g_hat_r[i] for i in range(half)]
    logger_cs_mpc.debug("Calculate Q_prime.")
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
        logger_cs_mpc.debug(f"Provers opened z_prime")
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
    logger_cs_mpc.debug(f"Provers opened y.")

    # Non-interactive proof
    # Step 1: Prover calculates ("announcement") t, A
    order = gf.order
    r = list(mpc._random(secfld) for i in range(n))
    rho = mpc._random(secfld)
    t = L(r)
    logger_cs_mpc.debug("Calculate A.")
    A = await vector_commitment(r, rho, g, h)
    t = await mpc.output(t)
    logger_cs_mpc.debug(f"Provers opened t, A.")
    proof["t"] = t
    proof["A"] = A

    # Step 2: Prover computes challenge
    # input_list = [t, A, generators, P, L, y]
    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.to_affine(), generators, P.to_affine(), L, y]
    else:
        input_list = [t, A, generators, P, L, y]

    logger_cs_mpc_hin.debug(
        f"Method protocol_5_prover: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c0 = pivot.fiat_shamir_hash(
        input_list + [0] + ["First hash of compressed pivot"], order
    )
    c1 = pivot.fiat_shamir_hash(
        input_list + [1] + ["First hash of compressed pivot"], order
    )
    logger_cs_mpc_hout.debug(f"After hash, hash=\n{c0}, {c1}")

    # Step 3: Prover calculates
    z = [c0 * x_i + r[i] for i, x_i in enumerate(x)]
    phi = c0 * gamma + rho
    z_hat = z + [phi]

    # Step 4: Prover calculates following public variables
    g_hat = g + [h]
    logger_cs_mpc.debug("Calculate Q.")
    Q = A * (P ** c0) * (k ** _int(c1 * (c0 * y + t)))
    L_tilde = pivot.LinearForm(L.coeffs + [0]) * c1
    # assert L(z)*c1 == L_tilde(z_hat)
    proof = await protocol_4_prover(g_hat, k, Q, L_tilde, z_hat, gf, proof)
    return proof


def calculate_fgh_polys(a, b, c, gf, secfld):
    # Calculate random polynomials f, g, h
    r_a = mpc._random(secfld)
    r_b = mpc._random(secfld)
    logger_cs_mpc.debug("Calculate f_poly.")
    f_poly = qc.Poly(qc.lagrange_interp_ff(a + [r_a], gf))
    logger_cs_mpc.debug("Calculate g_poly.")
    g_poly = qc.Poly(qc.lagrange_interp_ff(b + [r_b], gf))
    logger_cs_mpc.debug("Calculate h_poly.")
    h_poly = f_poly * g_poly
    logger_cs_mpc.debug("Done calculating f, g, h.")
    # assert c == [h_poly.eval(i+1) for i in range(m)], "Evaluations of h at 1..m not equal to vector c"
    return f_poly, g_poly, h_poly


async def protocol_8_excl_pivot_prover(generators, code, x, gf, use_koe=False):
    """ Non-interactive implementation of Protocol 8,
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

    logger_cs_mpc.debug(f"Calculate witness.")
    n = len(x)
    xc = cs.calculate_witness(code, x)
    proof = {}
    logger_cs_mpc.debug(f"Transform code to flatcode (circuit).")
    flatcode, inputs, varnames, r1cs = cs.code_to_flatcode_and_r1cs(code)

    mul_indices_of_flatcode = cs.mul_in_flatcode(flatcode)
    m = len(mul_indices_of_flatcode)
    output_variables = [symbol for symbol in varnames if symbol.startswith("~out")]

    # Calculate a, b, c vectors, only for mul gates
    logger_cs_mpc.debug(f"Calculate a, b, c vectors.")
    a, b = cs.calculate_ab_vectors(r1cs, xc, mul_indices_of_flatcode)
    c = mpc.schur_prod(a, b)

    # # Calculate random polynomials f, g, h
    # f_poly, g_poly, h_poly = calculate_fgh_polys(a, b, c, gf, secfld)

    # # Prover commits to (x, f(0), g(0), h(0), h(1), .., h(2m))
    # h_evaluations = [h_poly.eval(i + 1) for i in range(2 * m)]
    # z = x + [f_poly.eval(0), g_poly.eval(0), h_poly.eval(0)] + h_evaluations

    logger_cs_mpc.debug(f"Calculate z.")
    f0 = mpc._random(secfld)
    g0 = mpc._random(secfld)
    a = [f0] + a
    b = [g0] + b
    logger_cs_mpc.debug(f"Calculate z: a, b = await mpc.gather(a,b)")
    a, b = await mpc.gather(a, b)
    logger_cs_mpc.debug(f"Calculate z: fs = recombine(...)")
    fs = recombine(gf, list(zip(range(m+1), a)), list(range(m+1, 2*m+1)))
    logger_cs_mpc.debug(f"Calculate z: gs = recombine(...)")
    gs = recombine(gf, list(zip(range(m+1), b)), list(range(m+1, 2*m+1)))
    logger_cs_mpc.debug(f"Calculate z: hs = list(...)")
    hs = list(map(secfld, await mpc.schur_prod(fs, gs)))
    z = x + [f0, g0, f0 * g0] + c + hs

    gamma = mpc._random(secfld)

    if use_koe:
        S = range(len(z))
        z_commitment_P, z_commitment_pi = await koe_restriction_argument_prover(
            S, z, gamma, pp
        )
        z_commitment = {"P": z_commitment_P, "pi": z_commitment_pi}
        proof["z_commitment"] = z_commitment
    else:
        logger_cs_mpc.debug("Calculate commitment for z, denoted by [z].")
        z_commitment = await vector_commitment(z, gamma, g, h)
        # logger_cs_mpc.debug("Calculated [z] in secret shared form.")
        logger_cs_mpc.debug(f"Opened [z] (z_commitment).")
        proof["z_commitment"] = z_commitment

    logger_cs_mpc.debug("Apply first hash of circuit satisfiability protocol.")
    input_list = [z_commitment, code, "First hash circuit satisfiability protocol"]
    logger_cs_mpc_hin.debug(
        f"Method protocol_8_excl_pivot_prover (1): input_list={input_list}"
    )
    c = pivot.fiat_shamir_hash(input_list, gf.order)     # TODO: name clash c = a * b from above
    logger_cs_mpc_hout.debug(f"After hash, hash=\n{c}")

    # Prover's response
    # y1 = recombine(gf, list(zip(range(m+1), b)), c)  # f_poly.eval(c)  # TODO: typo!?? b -> a?
    # y2 = recombine(gf, list(zip(range(m+1), b)), c)  # g_poly.eval(c)
    # y3 = y1 * y2
    # NB: these values are not used below! TS: These previous lines (y1, y2, y3) are for debugging/checking consistency only.

    # Create linear forms corresponding to f(c), g(c), h(c) for Pi_Nullity step in Protocol 6 (Protocol 8 in updated version)
    linform_f, linform_g, linform_h = cs.create_fgh_linear_forms(
        r1cs, c, varnames, flatcode, mul_indices_of_flatcode, n, m, gf
    )

    y1 = linform_f(z)
    y2 = linform_g(z)
    y3 = linform_h(z)
    y1 = await mpc.output(y1)
    y2 = await mpc.output(y2)
    y3 = await mpc.output(y3)
    assert y1 * y2 == y3
    proof["y1"] = y1
    proof["y2"] = y2
    proof["y3"] = y3

    circuits = []
    outputs = []
    for output_var in output_variables:
        circuit = cs.express_as_x_or_gamma(output_var, flatcode, varnames, n)
        y = circuit(z)
        y = await mpc.output(y)
        # assert y == xc[varnames.index(output_var)], f"Output of circuit {y} not equal to ~out in witness."
        circuits.append(circuit)
        outputs.append(y)
    proof["outputs"] = outputs

    lin_forms = [circuit - y for circuit, y in zip(circuits, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    # Follow prove_nullity_compressed but add more inputs to Fiat-Shamir hash
    logger_cs_mpc.debug("Apply second hash of circuit satisfiability protocol.")
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
    logger_cs_mpc_hin.debug(
        f"Method protocol_8_excl_pivot_prover (2): input_list={input_list}"
    )
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_mpc_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    proof["L"] = L
    return proof, z_commitment, L, z, gamma


async def prove_linear_form_eval(g, h, P, L, y, x, gamma, gf):
    """ Sigma protocol Pi_s (protocol 2) from AC20.

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
    logger_cs_mpc.debug(f"Provers opened t.")
    logger_cs_mpc.debug(f"Provers opened A.")

    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.to_affine(), g, h, P.to_affine(), L, y]
    else:
        input_list = [t, A, g, h, P, L, y]
    logger_cs_mpc_hin.debug(f"Method prove_linear_form_eval: input_list={input_list}.")
    c = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_mpc_hout.debug(f"After hash, hash=\n{c}")
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
    """ Non-interactive implementation of Protocol 8, prover-side,
    including Nullity using compressed pivot (Protocol 5).
    """
    logger_cs_mpc.debug(f"Enter circuit_sat_prover. pivot_choice={pivot_choice}")

    logger_cs_mpc.debug(f"Start protocol 8, excluding pivot proof.")
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
        pivot_proof, u = await koe_opening_linear_form_prover(L, z, gamma, generators, P, pi)
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
