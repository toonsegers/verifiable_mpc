"""Demo Pinocchio protocol based on pynocchio module.

Credits:
* Trinocchio by Berry Schoenmakers, Meilof Veeningen 
  and Niels de Vreede: https://eprint.iacr.org/2015/480
* MPyC by Berry Schoenmakers: https://github.com/lschoe/mpyc
* Original bn256 module by Jack Lloyd: 
  https://github.com/randombit/pairings.py/ (BSD-2-Clause license)
* r1cs and qap modules by Vitalik Buterin: 
  https://github.com/ethereum/research/tree/master/zksnark (MIT license)
"""

import logging
import pprint
from random import SystemRandom
import os, sys

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import verifiable_mpc.trinocchio.pynocchio as pynocchio
from mpyc.finfields import GF, FiniteFieldElement
from mpyc.fingroups import EllipticCurve
import verifiable_mpc.tools.code_to_qap as c2q
import verifiable_mpc.tools.qap_creator as qc


def main():
    pp = pprint.PrettyPrinter(indent=4)

    bn_curve = EllipticCurve('BN256', 'jacobian')
    g1 = bn_curve.generator
    bn_twist = EllipticCurve('BN256_twist', 'jacobian')
    g2 = bn_twist.generator
    bn_curve.is_additive = True  
    bn_curve.is_multiplicative = False
    bn_twist.is_additive = True
    bn_twist.is_multiplicative = False
    # TODO: reconsider global setting of is_additive
    # interaction with other demos, e.g., in python -m unittest
    # necessitates resetting to default settings at this point

    modulus = bn_curve.order

    gf = GF(modulus=modulus)
    gf.is_signed = False

    # Set inputs
    inputs = [gf(3)]
    code = """
def qeval(x):
    y = x**3 + x**2 + x
    return y + x + 5
"""

    #     inputs = [gf(3), gf(2)]
    #     code = """
    # def qeval(x, y):
    #     z = x**3 + 2*y**2
    #     return z + x + 5
    # """

    # QAP creation step
    qap = c2q.QAP(code, gf)
    print(f"QAP created. Size: {qap.m}, degree {qap.d}.")

    # Trusted Party's KeyGen step
    td = pynocchio.Trapdoor(modulus)
    gen = pynocchio.Generators(td, g1, g2)
    evalkey = pynocchio.generate_evalkey(td, qap, gen)
    verikey = pynocchio.generate_verikey(td, qap, gen)
    print("Trusted setup completed.")

    # Prover's steps
    c = qap.calculate_witness(inputs)
    p = pynocchio.compute_p_poly(qap, c)
    h, r = p / qap.t
    assert r == qc.Poly(
        [0] * qap.d
    ), "Remainder of p(x)/t(x) for given witness is not 0"
    deltas = pynocchio.SampleDeltas(modulus)
    # Create zero knowledge variant of h poly
    h_zk = h + pynocchio.compute_h_zk_terms(qap, c, deltas)
    h = h_zk
    # Compute proof
    proof = pynocchio.compute_proof(qap, c, h, evalkey, deltas)
    print("Proof computed.")

    # Verifier's step
    verifications = pynocchio.verify(qap, verikey, proof, c[: qap.out_ix + 1])
    if all(verifications.values()):
        print("All checks passed.")
    else:
        print("Not all checks passed.")
    pp.pprint(verifications)
    return verifications    


if __name__ == "__main__":
    verifications = main()
#     # Globals when running script
#     pp = pprint.PrettyPrinter(indent=4)

#     bn_curve = EllipticCurve(ell.BN256, ell.WEI_JAC, ell.Weierstr_Jacobian_Arithm)
#     g1 = bn_curve.base_pt
#     bn_twist = EllipticCurve(ell.BN256_TWIST, ell.WEI_JAC, ell.Weierstr_Jacobian_Arithm)
#     g2 = bn_twist.base_pt

#     modulus = bn_curve.order

#     gf = GF(modulus=modulus)
#     gf.is_signed = False

#     # Set inputs
#     inputs = [gf(3)]
#     code = """
# def qeval(x):
#     y = x**3 + x**2 + x
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
#     td = pynocchio.Trapdoor(modulus)
#     gen = pynocchio.Generators(td, g1, g2)
#     evalkey = pynocchio.generate_evalkey(td, qap, gen)
#     verikey = pynocchio.generate_verikey(td, qap, gen)
#     print("Trusted setup completed.")

#     # Prover's steps
#     c = qap.calculate_witness(inputs)
#     p = pynocchio.compute_p_poly(qap, c)
#     h, r = p / qap.t
#     assert r == qc.Poly(
#         [0] * qap.d
#     ), "Remainder of p(x)/t(x) for given witness is not 0"
#     deltas = pynocchio.SampleDeltas(modulus)
#     # Create zero knowledge variant of h poly
#     h_zk = h + pynocchio.compute_h_zk_terms(qap, c, deltas)
#     h = h_zk
#     # Compute proof
#     proof = pynocchio.compute_proof(qap, c, h, evalkey, deltas)
#     print("Proof computed.")

#     # Verifier's step
#     verifications = pynocchio.verify(qap, verikey, proof, c[: qap.out_ix + 1])
#     if all(verifications.values()):
#         print("All checks passed.")
#     else:
#         print("Not all checks passed.")
#     pp.pprint(verifications)


# TODO list
# TODO-1: Replace simple classes by SimpleNamespaces: https://stackoverflow.com/questions/16279212/how-to-use-dot-notation-for-dict-in-python
# TODO-2: Optimize g_eval: consider caching, or reuse of v_g1 set for alpha v_g1 set
# TODO-3: Pinocchio Section 4.2.1 (p.8) discusses bad asymptotics with naive algorithms instead of FFT; update with FFT
