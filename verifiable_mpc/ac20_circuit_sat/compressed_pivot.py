""" Implementation of https://eprint.iacr.org/2020/152
``Compressed Σ-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

Protocols:
* Compressed Σ-protocol Π_c for relation R, page 15, protocol 5
"""

import os
import sys
import logging
from random import SystemRandom

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import verifiable_mpc.ac20_circuit_sat.pivot as pivot
from sec_groups.fingroups import EllipticCurveElement


prng = SystemRandom()

logger_cp = logging.getLogger("compressed_pivot")
logger_cp.setLevel(logging.INFO)

logger_cp_hin = logging.getLogger("compressed_pivot_hash_inputs")
logger_cp_hin.setLevel(logging.INFO)

logger_cp_hout = logging.getLogger("compressed_pivot_hash_outputs")
logger_cp_hout.setLevel(logging.INFO)


def protocol_4_prover(g_hat, k, Q, L_tilde, z_hat, gf, proof={}, round_i=0):
    """ Non-interactive version of Protocol 4 from Section 4.2: Prover's side

    """

    # Step 5: Prover calculates A, B
    half = len(g_hat) // 2
    g_hat_l = g_hat[:half]
    g_hat_r = g_hat[half:]
    z_hat_l = z_hat[:half]
    z_hat_r = z_hat[half:]
    logger_cp.debug("Calculate A_i, B_i.")
    A = pivot.vector_commitment(z_hat_l, int(L_tilde([0] * half + z_hat_l)), g_hat_r, k)
    B = pivot.vector_commitment(z_hat_r, int(L_tilde(z_hat_r + [0] * half)), g_hat_l, k)

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

    logger_cp_hin.debug(
        f"Method protocol_4_prover: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c = pivot.fiat_shamir_hash(input_list, order)
    logger_cp_hout.debug(f"After hash, hash=\n{c}")

    # Step 7: Prover calculates following public values
    logger_cp.debug("Calculate g_prime.")
    g_prime = [(g_hat_l[i] ** c) * g_hat_r[i] for i in range(half)]
    logger_cp.debug("Calculate Q_prime.")
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
        proof["z_prime"] = z_prime
        return proof
    else:
        # Step 9: Prover sends z_prime and starts recursion of protocol 4
        round_i += 1
        proof = protocol_4_prover(
            g_prime, k, Q_prime, L_prime, z_prime, gf, proof, round_i
        )
        return proof


def protocol_5_prover(generators, P, L, y, x, gamma, gf):
    g = generators["g"]
    h = generators["h"]
    k = generators["k"]

    proof = {}
    # Public parameters
    n = len(x)
    L, y = pivot.affine_to_linear(L, y, n)
    assert (
        bin(n + 1).count("1") == 1
    ), "This implementation requires n+1 to be power of 2 (else, use padding with zeros)."

    # Non-interactive proof
    # Step 1: Prover calculates ("announcement") t, A
    order = gf.order
    r = list(prng.randrange(order) for i in range(n))
    rho = prng.randrange(order)
    logger_cp.debug("Calculate t.")
    t = L(r)
    logger_cp.debug("Calculate A.")
    A = pivot.vector_commitment(r, rho, g, h)

    proof["t"] = t
    proof["A"] = A

    # Step 2: Prover computes challenge
    # input_list = [t, A, generators, P, L, y]
    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.to_affine(), generators, P.to_affine(), L, y]
    else:
        input_list = [t, A, generators, P, L, y]

    logger_cp_hin.debug(
        f"Method protocol_5_prover: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c0 = pivot.fiat_shamir_hash(
        input_list + [0] + ["First hash of compressed pivot"], order
    )
    c1 = pivot.fiat_shamir_hash(
        input_list + [1] + ["First hash of compressed pivot"], order
    )
    logger_cp_hout.debug(f"After hash, hash=\n{c0}, {c1}")

    # Step 3: Prover calculates
    z = [c0 * x_i + r[i] for i, x_i in enumerate(x)]
    phi = gf(c0 * gamma + rho)
    z_hat = z + [phi]
    # Step 4: Prover calculates following public variables
    g_hat = g + [h]
    logger_cp.debug("Calculate Q.")
    Q = A * (P ** c0) * (k ** int(c1 * (c0 * y + t)))
    L_tilde = pivot.LinearForm(L.coeffs + [0]) * c1
    assert L(z) * c1 == L_tilde(z_hat)

    proof = protocol_4_prover(g_hat, k, Q, L_tilde, z_hat, gf, proof)
    return proof


def protocol_4_verifier(g_hat, k, Q, L_tilde, gf, proof, round_i=0):
    """ Non-interactive version of Protocol 4 from Section 4.2: Verifier's side

    """

    # Step 5
    half = len(g_hat) // 2
    g_hat_l = g_hat[:half]
    g_hat_r = g_hat[half:]

    logger_cp.debug("Load from proof: A_i, B_i.")
    A = proof["A" + str(round_i)]
    B = proof["B" + str(round_i)]

    # Step 6
    order = k.order
    # Hash A, B and all public elements of Protocol 4
    # input_list = [A, B, g_hat, k, Q, L_tilde]
    if isinstance(A, EllipticCurveElement):
        input_list = [A.to_affine(), B.to_affine(), g_hat, k, Q.to_affine(), L_tilde]
    else:
        input_list = [A, B, g_hat, k, Q, L_tilde]
    logger_cp_hin.debug(
        f"Method protocol_4_verifier: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c = pivot.fiat_shamir_hash(input_list, order)
    logger_cp_hout.debug(f"After hash, hash=\n{c}")

    # Step 7
    logger_cp.debug("Calculate g_prime.")
    g_prime = [(g_hat_l[i] ** c) * g_hat_r[i] for i in range(half)]
    logger_cp.debug("Calculate Q_prime.")
    Q_prime = A * (Q ** c) * (B ** (c ** 2))

    assert (
        L_tilde.constant == 0
    ), "Next line assumes L_tilde is a linear form, not affine form."
    c_L_tilde_l_coeffs = [coeff * gf(c) for coeff in L_tilde.coeffs[:half]]
    L_prime = pivot.LinearForm(c_L_tilde_l_coeffs) + pivot.LinearForm(
        L_tilde.coeffs[half:]
    )

    # Step 8
    if len(g_prime) <= 2:
        z_prime = proof["z_prime"]
        Q_check = pivot.vector_commitment(z_prime, int(L_prime(z_prime)), g_prime, k)
        logger_cp.debug("Arrived in final step of protocol_4_verifier.")
        logger_cp.debug(f"Q_check= {Q_check}")
        logger_cp.debug(f"Q_prime= {Q_prime}")
        verification = Q_check == Q_prime
        return verification
    else:
        # Step 9
        round_i += 1
        return protocol_4_verifier(g_prime, k, Q_prime, L_prime, gf, proof, round_i)


def protocol_5_verifier(generators, P, L, y, proof, gf):
    g = generators["g"]
    h = generators["h"]
    k = generators["k"]

    order = gf.order
    n = len(g)
    L, y = pivot.affine_to_linear(L, y, n)
    logger_cp.debug("Load from proof: t, A.")
    t = proof["t"]
    A = proof["A"]
    # input_list = [t, A, generators, P, L, y]
    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.to_affine(), generators, P.to_affine(), L, y]
    else:
        input_list = [t, A, generators, P, L, y]

    logger_cp_hin.debug(
        f"Method protocol_5_verifier: Before fiat_shamir_hash, input_list=\n{input_list}"
    )
    c0 = pivot.fiat_shamir_hash(
        input_list + [0] + ["First hash of compressed pivot"], order
    )
    c1 = pivot.fiat_shamir_hash(
        input_list + [1] + ["First hash of compressed pivot"], order
    )
    logger_cp_hout.debug(f"After hash, hash=\n{c0}, {c1}")

    g_hat = g + [h]
    logger_cp.debug("Calculate Q.")
    Q = A * (P ** c0) * (k ** int(c1 * (c0 * y + t)))
    L_tilde = pivot.LinearForm(L.coeffs + [0]) * c1

    verification = protocol_4_verifier(g_hat, k, Q, L_tilde, gf, proof)
    return verification


# TODO-list
# TODO1: Input y can be removed from inputs to protocol5 method.
