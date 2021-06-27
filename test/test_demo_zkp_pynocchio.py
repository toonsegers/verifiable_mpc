""" Tests for pynocchio.py 

"""

import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import unittest

from demos import demo_zkp_pynocchio


class CircuitSat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_pynocchio(self):
        verification = demo_zkp_pynocchio.main()
        self.assertEqual(all(verification.values()), True)


if __name__ == "__main__":
    unittest.main()


# import pprint
# from random import SystemRandom

# import bn256
# import code_to_qap as c2q
# import trinocchi.pynocchio
# import tools.code_to_r1cs_and_qap.qap_creator as qc
# from mpyc.finfields import GF, FiniteFieldElement

# pp = pprint.PrettyPrinter(indent=4)

# code = """
# def qeval(x):
#     y = x**3
#     return y + x + 5
# """

# if __name__ == '__main__':
#     gf = GF(modulus = bn256.order)
#     gf.is_signed = False
#     prng = SystemRandom()
#     g1 = bn256.curve_G
#     g2 = bn256.twist_G

#     print("Test polynomial interpolation over finite field (randomized)...")
#     degree = 10
#     poly = [gf(prng.randrange(bn256.order)) for i in range(degree+1)]
#     y = [qc.eval_poly(poly, i+1) for i in range(degree+1)]
#     poly2 = qc.lagrange_interp_ff(y, gf)
#     assert poly == poly2, "Polynomial interpolation error"

#     print("Test correctness of target polynomial (randomized)...")
#     r1 = prng.randint(0, bn256.order)
#     r2 = prng.randint(0, bn256.order)
#     r3 = prng.randint(0, bn256.order)
#     code = f"""
# def qeval(x):
#     y = x**3 + {gf(r1)}*x**2 + {gf(r2)}*x
#     return y + x + {gf(r3)}
# """
#     qap = c2q.QAP(code, gf)
#     y = [qc.eval_poly(qap.t, i) for i in range(1, qap.d+1)]
#     assert y == [0]*qap.d, "Target polynomial is not zero for all roots corresponding to gates"

#     print("Test if trapdoor contains None...")
#     td = pynocchio.Trapdoor(bn256.order)
#     assert(None not in [td.r_v, td.r_w, td.r_y, td.s, td.alpha_v, td.alpha_w, td.alpha_y, td.beta, td.gamma]), "Trapdoor contains None"

#     print("Test g_eval method...")
#     assert str(g2.scalar_mul(0)) == str(pynocchio.g_eval(g2, qap.t, 4, td.alpha_v)), "g^(alpha*t(4)) != id; should be identity element"
#     vs = qc.eval_poly(qap.v[3], td.s)
#     g1vs = g1.scalar_mul(int(vs))
#     assert str(pynocchio.g_eval(g1, qap.v[3], td.s)) == str(g1vs)

#     print("Test qap indices...")
#     print(f"{qap.varnames=}")
#     print(f"{len(qap.varnames)=}")
#     print(f"{qap.m=}")
#     print(f"{qap.d=}")
#     print(f"{qap.indices_io=}")
#     print(f"{qap.indices_mid=}")
#     print(f"{qap.indices=}")
