""" Implementation of https://eprint.iacr.org/2020/152
``Compressed Σ-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

Constant-sized, non-interactive proof based on the n-power
Knowledge-of-Exponent Assumption (n-PKEA).

Protocols:
* Program from the Knowledge-of-Exponent Assumption,
Section 9, page 25
"""

import os
import sys
from random import SystemRandom

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import mpyc.mpctools as mpctools
import verifiable_mpc.ac20_circuit_sat.pivot as pivot
import verifiable_mpc.tools.qap_creator as qc
import sec_groups.pairing as pairing

prng = SystemRandom()



def list_mul(x):
    return mpctools.reduce(type(x[0]).operation, x)


def vector_commitment(x, gamma, g, h):
    """ Pedersen vector commitment. Definition 1 of AC20.
    """
    assert len(g) >= len(x), "Not enough generators."
    gpowers = [g[i] ** int(x_i) for i, x_i in enumerate(x)]
    # prod = mpctools.reduce(type(h).operation, gpowers)
    prod = list_mul(gpowers)
    
    c = (h ** gamma) * prod
    return c


def _pairing(a, b):
    """Flip inputs to pairing.optimal_ate() function.

    First input "base" curve, second input "twist" curve 
    as per Pinocchio/Trinocchio notation.
    """
    return pairing.optimal_ate(b, a)


def trusted_setup(_g1, _g2, n, order, progress_bar=False):
    g_exp = prng.randrange(1, order)
    alpha = prng.randrange(order)
    z = prng.randrange(order)
    g1 = _g1 ** g_exp
    g2 = (_g2 ** g_exp) ** alpha

    pp_lhs = []
    pp_rhs = []
    g1_base = g1
    g2_base = g2
    for i in range(2 * n):
        g1_base = g1_base ** z
        g2_base = g2_base ** z

        pp_lhs.append(g1_base)
        pp_rhs.append(g2_base)

        if progress_bar:
            print(f"Generating keys: {round(100*i/(2*n+1))}%", end="\r")

    pp = {"pp_lhs": pp_lhs, "pp_rhs": pp_rhs}
    return pp


def restriction_argument_prover(S, x, gamma, pp):
    """Resriction argument from [Gro10]: Prover's side.

    Extension of the Knowledge Commitment Scheme, but proofs that a subset
    of indices S of {0, ..., n-1}, corresponding to vector x, was used.
    """

    # Prover commits only to S-indices of x with (new) commitment P, with secret random gamma
    P = (pp["pp_lhs"][0] ** gamma) * list_mul([pp["pp_lhs"][i + 1] ** x[i] for i in S])

    # Trusted party or verifier samples beta; only required in designated verifier scenario
    # beta = prng.randrange(1, order)
    # sigma = (g1**beta, [(g2**(z**i))**beta for i in S])

    # Prover computes pi^alpha (slight deviation from Thomas' notes, who computes pi)
    pi = (pp["pp_rhs"][0] ** gamma) * list_mul([pp["pp_rhs"][i + 1] ** x[i] for i in S])
    return P, pi


def restriction_argument_verifier(P, pi, pp):
    """ Resriction argument from [Gro10]: Verifier's side
    """
    verification = _pairing(P, pp["pp_rhs"][0]) == _pairing(pp["pp_lhs"][0], pi)
    return verification


def opening_linear_form_prover(L, x, gamma, pp, P=None, pi=None):
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
        P, pi = restriction_argument_prover(S, x, gamma, pp)
    proof["P"] = P
    proof["pi"] = pi

    # Prover computes coefficients of c_poly and Q, sends Q to verifier
    u = L(x)
    L_linear, u_linear = pivot.affine_to_linear(L, u, n)

    c_poly_lhs = qc.Poly([gamma] + [x_i for x_i in x])
    c_poly_rhs = qc.Poly([L_linear.coeffs[n - (j + 1)] for j in range(n)])
    c_poly = c_poly_lhs * c_poly_rhs

    assert u_linear == c_poly.coeffs[n], "L(x) not equal to n-th coefficient of c_poly"
    c_bar = c_poly.coeffs
    c_bar[n] = 0
    assert len(pp["pp_lhs"]) == 2 * n
    Q = list_mul([g_i ** (-1 * c_bar[i]) for i, g_i in enumerate(pp["pp_lhs"])])
    proof["Q"] = Q
    return proof, u


def opening_linear_form_verifier(L, pp, proof, u):
    n = len(L.coeffs)
    g1 = pp["pp_lhs"][0]
    g2 = pp["pp_rhs"][0]
    L_linear, u_linear = pivot.affine_to_linear(L, u, n)
    P = proof["P"]
    pi = proof["pi"]
    Q = proof["Q"]
    verification = {}
    verification["restriction_arg_check"] = _pairing(P, g2) == _pairing(g1, pi)
    R = list_mul([pp["pp_rhs"][j] ** (L_linear.coeffs[n - (j + 1)]) for j in range(n)])
    check_lhs = _pairing(P, R) * _pairing(Q, g2)
    check_rhs = _pairing(g1, pp["pp_rhs"][n] ** u_linear)
    verification["PRQ_check"] = check_lhs == check_rhs
    return verification


def prove_nullity_koe(pp, lin_forms, x, gamma, gf, P, pi):
    """ Nullity protocol calling the ZK Proof for lienar form openings
    based on Knowledge of Exponent assumption.
    """
    input_list = [P, lin_forms]
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    L = pivot.LinearForm([gf(c) if isinstance(c, int) else c for c in L.coeffs])
    proof, u = opening_linear_form_prover(L, x, gamma, pp, P, pi)
    return proof, L, u
