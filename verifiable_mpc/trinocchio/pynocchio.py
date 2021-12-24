"""Implementation of Pinocchio [PGHR13] Protocol 2.

Paper by Parno, Gentry, Howell, Raykova (2013):
https://eprint.iacr.org/2013/279
Python implementation using BN256 curve and asymmetric bilinear map.

Small adaptations to original Pinocchio protocol:
* H-check and construction of h-polynomial in zero knowledge case
  follows Trinocchio: https://eprint.iacr.org/2015/480
* Compute s-powers in evalkey with range(0, qap.d + 1) instead of
  range(1, qap.d + 1), same for h_g1_terms in compute_proof()

Credits:
* Original bn256 module by Jack Lloyd:
  https://github.com/randombit/pairings.py/ (BSD-2-Clause license)
* Original r1cs and qap tools by Vitalik Buterin:
  https://github.com/ethereum/research/tree/master/zksnark (MIT license)
"""

import os, sys
project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import pprint
from random import SystemRandom

from mpyc.finfields import GF
import verifiable_mpc.tools.code_to_qap as c2q
import verifiable_mpc.tools.qap_creator as qc
from verifiable_mpc.ac20.pairing import optimal_ate
from mpyc.fingroups import FiniteGroupElement


# Globals for module
prng = SystemRandom()
# Point point_add to generic group operation
point_add = FiniteGroupElement.__matmul__

class Trapdoor:
    def __init__(self, modulus):
        _td = list(prng.randrange(modulus) for i in range(8))
        r_v, r_w, s, alpha_v, alpha_w, alpha_y, beta, gamma = _td
        r_y = r_v * r_w % modulus
        self.r_v = r_v
        self.r_w = r_w
        self.r_y = r_y
        self.s = s
        self.alpha_v = alpha_v
        self.alpha_w = alpha_w
        self.alpha_y = alpha_y
        self.beta = beta
        self.gamma = gamma


class SampleDeltas:
    def __init__(self, modulus):
        _deltas = list(prng.randrange(modulus) for i in range(3))
        delta_v, delta_w, delta_y = _deltas
        self.v = delta_v
        self.w = delta_w
        self.y = delta_y


class Generators:
    def __init__(self, td, g1, g2):
        self.g1 = g1
        self.g2 = g2
        self.g1_v = td.r_v * g1
        self.g1_w = td.r_w * g1
        self.g2_w = td.r_w * g2
        self.g1_y = td.r_y * g1
        self.g2_y = td.r_y * g2


def pairing(a, b):
    """Perform pairing operation.

    First input "base" curve, second input "twist" curve as per Pinocchio/Trinocchio notation.
    """
    return optimal_ate(b, a)




def apply_to_list(op, inputs):
    """Apply operator `op` to list `inputs` via binary tree
    """
    n = len(inputs)
    if n == 1:
        return inputs[0]

    m0 = apply_to_list(op, inputs[: n // 2])
    m1 = apply_to_list(op, inputs[n // 2 :])
    return op(m0, m1)


def g_eval(gen, poly, s, alpha=1):
    """ Evaluate poly at s, then scale curve point with outcome
    """
    poly_at_s = poly.eval(s)
    return int(alpha * poly_at_s) * gen


def generate_evalkey(td, qap, gen):
    """ Generate public evaluation key.

    Note that g_w^w(s) terms are on twist curve (with generator g2) (See Pinocchio, p. 6)
    """
    v_g1 = {
        "r_v*v" + str(i) + "*g1": g_eval(gen.g1_v, qap.v[i], td.s)
        for i in qap.indices_mid
    }
    w_g2 = {
        "r_w*w" + str(i) + "*g2": g_eval(gen.g2_w, qap.w[i], td.s)
        for i in qap.indices_mid
    }
    y_g1 = {
        "r_y*y" + str(i) + "*g1": g_eval(gen.g1_y, qap.y[i], td.s)
        for i in qap.indices_mid
    }
    alphav_g1 = {
        "r_v*alpha_v*v" + str(i) + "*g1": g_eval(gen.g1_v, qap.v[i], td.s, td.alpha_v)
        for i in qap.indices_mid
    }
    alphaw_g1 = {
        "r_w*alpha_w*w" + str(i) + "*g1": g_eval(gen.g1_w, qap.w[i], td.s, td.alpha_w)
        for i in qap.indices_mid
    }
    alphay_g1 = {
        "r_y*alpha_y*y" + str(i) + "*g1": g_eval(gen.g1_y, qap.y[i], td.s, td.alpha_y)
        for i in qap.indices_mid
    }
    spowers_g1 = {
        "s^" + str(i) + "*g1": td.s**i * gen.g1 for i in range(0, qap.d + 1)
    }
    betavwy_g1 = {
        "r_v*beta*v+r_w*beta*w+r_y*beta*y"
        + str(i)
        + "_g1": g_eval(gen.g1_v, qap.v[i], td.s, td.beta)
        + (g_eval(gen.g1_w, qap.w[i], td.s, td.beta))
        + (g_eval(gen.g1_y, qap.y[i], td.s, td.beta))
        for i in qap.indices_mid
    }

    # Elements required to make proofs zero knowledge
    zk_elements = {
        "r_v*t*g1": g_eval(gen.g1_v, qap.t, td.s),
        "r_w*t*g2": g_eval(gen.g2_w, qap.t, td.s),
        "r_y*t*g1": g_eval(gen.g1_y, qap.t, td.s),
        "r_v*alpha_v*t*g1": g_eval(gen.g1_v, qap.t, td.s, td.alpha_v),
        "r_w*alpha_w*t*g1": g_eval(gen.g1_w, qap.t, td.s, td.alpha_w),
        "r_y*alpha_y*t*g1": g_eval(gen.g1_y, qap.t, td.s, td.alpha_y),
        "r_v*beta*t*g1": g_eval(gen.g1_v, qap.t, td.s, td.beta),
        "r_w*beta*t*g1": g_eval(gen.g1_w, qap.t, td.s, td.beta),
        "r_y*beta*t*g1": g_eval(gen.g1_y, qap.t, td.s, td.beta),
        "t*g1": g_eval(gen.g1, qap.t, td.s),
    }

    evalkey = {
        **v_g1,
        **w_g2,
        **y_g1,
        **alphav_g1,
        **alphaw_g1,
        **alphay_g1,
        **spowers_g1,
        **betavwy_g1,
        **zk_elements,
    }
    return evalkey


def generate_verikey(td, qap, gen):
    """ Generate public verification key.

    """
    part1 = {
        "g1": gen.g1,
        "g2": gen.g2,
        "alpha_v*g2": td.alpha_v * gen.g2,
        "alpha_w*g1": td.alpha_w * gen.g1,
        "alpha_y*g2": td.alpha_y * gen.g2,
        "gamma*g2": td.gamma * gen.g2,
        "beta*gamma*g1": (td.beta * td.gamma) * gen.g1,
        "beta*gamma*g2": (td.beta * td.gamma) * gen.g2,
        "r_y*t*g2": g_eval(
            gen.g2_y, qap.t, td.s
        ),  # TS: added when following Trinocchio verification steps
    }
    v_g1 = {
        "r_v*v" + str(i) + "*g1": g_eval(gen.g1_v, qap.v[i], td.s)
        for i in qap.indices_io_and_0
    }
    w_g2 = {
        "r_w*w" + str(i) + "*g2": g_eval(gen.g2_w, qap.w[i], td.s)
        for i in qap.indices_io_and_0
    }
    y_g1 = {
        "r_y*y" + str(i) + "*g1": g_eval(gen.g1_y, qap.y[i], td.s)
        for i in qap.indices_io_and_0
    }
    verikey = {**part1, **v_g1, **w_g2, **y_g1}
    return verikey


def compute_p_poly(qap, c):
    """ Create QAP polynomial p, as per Definition 2 of Pinocchio, p. 3

    """
    v_terms = apply_to_list(qc.add_polys, [qap.v[i] * c[i] for i in qap.indices])
    w_terms = apply_to_list(qc.add_polys, [qap.w[i] * c[i] for i in qap.indices])
    y_terms = apply_to_list(qc.add_polys, [qap.y[i] * c[i] for i in qap.indices])
    p = (v_terms * w_terms) - y_terms
    return p


def compute_h_zk_terms(qap, c, deltas):
    """ Create QAP polynomial p, as per Definition 2 of Pinocchio, p. 3

    """
    v_terms = apply_to_list(
        qc.add_polys, [qap.w[i] * (c[i] * deltas.v) for i in qap.indices]
    )
    w_terms = apply_to_list(
        qc.add_polys, [qap.v[i] * (c[i] * deltas.w) for i in qap.indices]
    )
    h_zk_terms = v_terms + w_terms + qap.t * (deltas.v * deltas.w) - qc.Poly([deltas.y])
    return h_zk_terms


def compute_proof(qap, c, h, evalkey, deltas=None):
    vmid_g1_terms = [int(c[i]) * evalkey["r_v*v" + str(i) + "*g1"] for i in qap.indices_mid]
    wmid_g2_terms = [int(c[i]) * evalkey["r_w*w" + str(i) + "*g2"] for i in qap.indices_mid]
    ymid_g1_terms = [int(c[i]) * evalkey["r_y*y" + str(i) + "*g1"] for i in qap.indices_mid]
    alphavmid_g1_terms = [int(c[i]) * evalkey[f'r_v*alpha_v*v{i}*g1'] for i in qap.indices_mid]
    alphawmid_g1_terms = [int(c[i]) * evalkey[f'r_w*alpha_w*w{i}*g1'] for i in qap.indices_mid]
    alphaymid_g1_terms = [int(c[i]) * evalkey[f'r_y*alpha_y*y{i}*g1'] for i in qap.indices_mid]
    betavwymid_g1_terms = [int(c[i]) * evalkey[f'r_v*beta*v+r_w*beta*w+r_y*beta*y{i}_g1'] for i in qap.indices_mid]
    # h_g1_first_term = gen.g1 * h.coeffs[0]
    h_g1_terms = [int(h.coeffs[i]) * evalkey["s^" + str(i) + "*g1"] for i in range(0, len(h))]
    # h_g1_terms.append(h_g1_first_term)
    vmid_g1 = apply_to_list(point_add, vmid_g1_terms)
    wmid_g2 = apply_to_list(point_add, wmid_g2_terms)
    ymid_g1 = apply_to_list(point_add, ymid_g1_terms)
    alphavmid_g1 = apply_to_list(point_add, alphavmid_g1_terms)
    alphawmid_g1 = apply_to_list(point_add, alphawmid_g1_terms)
    alphaymid_g1 = apply_to_list(point_add, alphaymid_g1_terms)
    betavwymid_g1 = apply_to_list(point_add, betavwymid_g1_terms)
    h_g1 = apply_to_list(point_add, h_g1_terms)

    # In zero knowledge case, add terms required to make proof zero knowledge
    if deltas != None:
        vmid_g1 = vmid_g1 + deltas.v * evalkey["r_v*t*g1"]
        wmid_g2 = wmid_g2 + deltas.w * evalkey["r_w*t*g2"]
        ymid_g1 = ymid_g1 + deltas.y * evalkey["r_y*t*g1"]
        alphavmid_g1 = alphavmid_g1 + deltas.v * evalkey["r_v*alpha_v*t*g1"]
        alphawmid_g1 = alphawmid_g1 + deltas.w * evalkey["r_w*alpha_w*t*g1"]
        alphaymid_g1 = alphaymid_g1 + deltas.y * evalkey["r_y*alpha_y*t*g1"]
        betavwymid_g1 = (
            betavwymid_g1
            + deltas.v * evalkey["r_v*beta*t*g1"]
            + deltas.w * evalkey["r_w*beta*t*g1"]
            + deltas.y * evalkey["r_y*beta*t*g1"]
        )

    proof = {
        "r_v*v_mid*g1": vmid_g1,
        "r_w*w_mid*g2": wmid_g2,
        "r_y*y_mid*g1": ymid_g1,
        "r_v*alpha_v*v_mid*g1": alphavmid_g1,
        "r_w*alpha_w*w_mid*g1": alphawmid_g1,
        "r_y*alpha_y*y_mid*g1": alphaymid_g1,
        "r_v*beta*v_mid+r_w*beta*w_mid+r_y*beta*y_mid*g1": betavwymid_g1,
        "h*g1": h_g1,
    }
    return proof


def verify(qap, verikey, proof, c):
    verification = {}

    # Divisibility check
    vio_g1_terms = [int(c[i]) * verikey["r_v*v" + str(i) + "*g1"] for i in qap.indices_io]
    wio_g2_terms = [int(c[i]) * verikey["r_w*w" + str(i) + "*g2"] for i in qap.indices_io]
    yio_g1_terms = [int(c[i]) * verikey["r_y*y" + str(i) + "*g1"] for i in qap.indices_io]
    vio_g1 = apply_to_list(point_add, vio_g1_terms)
    wio_g2 = apply_to_list(point_add, wio_g2_terms)
    yio_g1 = apply_to_list(point_add, yio_g1_terms)
    lhs1 = pairing(
        verikey["r_v*v0*g1"] + vio_g1 + proof["r_v*v_mid*g1"],
        verikey["r_w*w0*g2"] + wio_g2 + (proof["r_w*w_mid*g2"]),
    )
    lhs2 = pairing(yio_g1 + proof["r_y*y_mid*g1"], verikey["g2"]).inverse()
    rhs = pairing(proof["h*g1"], verikey["r_y*t*g2"])
    lhs = lhs1 * lhs2
    verification["H"] = lhs == rhs

    # Linear combination check
    lhs = pairing(
        proof["r_v*v_mid*g1"], verikey["alpha_v*g2"]
    )  # following order as per Trinocchio p. 19
    rhs = pairing(proof["r_v*alpha_v*v_mid*g1"], verikey["g2"])
    verification["V"] = lhs == rhs

    lhs = pairing(
        verikey["alpha_w*g1"], proof["r_w*w_mid*g2"]
    )  # following order as per Trinocchio p. 19
    rhs = pairing(proof["r_w*alpha_w*w_mid*g1"], verikey["g2"])
    alphawmid_g1 = proof["r_w*alpha_w*w_mid*g1"]
    verikey_g2 = verikey["g2"]
    verification["W"] = lhs == rhs

    lhs = pairing(proof["r_y*alpha_y*y_mid*g1"], verikey["g2"])
    rhs = pairing(proof["r_y*y_mid*g1"], verikey["alpha_y*g2"])
    verification["Y"] = lhs == rhs

    # Check if same witness was used
    lhs = pairing(
        proof["r_v*beta*v_mid+r_w*beta*w_mid+r_y*beta*y_mid*g1"], verikey["gamma*g2"]
    )
    rhs1 = pairing(
        proof["r_v*v_mid*g1"] + (proof["r_y*y_mid*g1"]), verikey["beta*gamma*g2"]
    )
    rhs2 = pairing(verikey["beta*gamma*g1"], proof["r_w*w_mid*g2"])
    rhs = rhs1 * rhs2
    verification["Z"] = lhs == rhs

    return verification


# if __name__ == "__main__":
#     # Globals when running script
#     pp = pprint.PrettyPrinter(indent=4)
#     modulus = bn256.order
#     g1 = bn256.curve_G
#     g2 = bn256.twist_G
#     gf = GF(modulus=modulus)
#     gf.is_signed = False

#     # Set inputs
#     inputs = [gf(3)]
#     code = """
# def qeval(x):
#     y = x**3
#     return y + x + 5
# """

#     #     inputs = [gf(3), gf(2)]
#     #     code = """
#     # def qeval(x, y):
#     #     z = x**3 + 2*y**2
#     #     return z + x + 5
#     # """

#     # QAP creation step
#     qap = c2q.QAP(code, gf)
#     print(f"QAP created. Size: {qap.m}, degree {qap.d}.")

#     # Trusted Party's KeyGen step
#     td = Trapdoor(modulus)
#     gen = Generators(td, g1, g2)
#     evalkey = generate_evalkey(td, qap, gen)
#     verikey = generate_verikey(td, qap, gen)
#     print("Trusted setup completed.")

#     # Prover's steps
#     c = qap.calculate_witness(inputs)
#     p = compute_p_poly(qap, c)
#     h, r = p / qap.t
#     assert r == qc.Poly(
#         [0] * qap.d
#     ), "Remainder of p(x)/t(x) for given witness is not 0"
#     deltas = SampleDeltas(modulus)
#     # Create zero knowledge variant of h poly
#     h_zk = h + compute_h_zk_terms(qap, c, deltas)
#     h = h_zk
#     print(len(h.coeffs))
#     # Compute proof
#     proof = compute_proof(qap, c, h, evalkey, deltas)
#     print("Proof computed.")

#     # Verifier's step
#     verifications = verify(qap, verikey, proof, c[: qap.out_ix + 1])
#     if all(check == True for check in verifications.values()):
#         print("All checks passed.")
#     else:
#         print("Not all checks passed.")
#     pp.pprint(verifications)


# TODO list
# TODO-1: Replace simple classes by SimpleNamespaces: https://stackoverflow.com/questions/16279212/how-to-use-dot-notation-for-dict-in-python
# TODO-2: Optimize g_eval: consider caching, or reuse of v_g1 set for alpha v_g1 set
# TODO-3: Pinochcio Section 4.2.1 (p.8) discusses bad asymptotics with naive algorithms instead of FFT; update with FFT
