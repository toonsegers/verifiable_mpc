# # WORK IN PROGRESS, fix bug that came up after migration to sec_grps2

# import sys, os
# project_root = sys.path.append(os.path.abspath('../..'))
# if project_root not in sys.path:
#     sys.path.insert(0, project_root)

# from random import SystemRandom
# import unittest

# import ac20_circuit_sat.pivot as pivot
# import ac20_circuit_sat.knowledge_of_exponent as koe
# from mpyc.finfields import GF
# # import zk_helpers.pairings.bn256_to_multiplicative as bn256
# from sec_grps2.fingroups import QuadraticResidue, EllipticCurve
# import sec_grps2.ellcurves as ell

# class KoE(unittest.TestCase):

#     @classmethod
#     def setUpClass(cls):
#         cls.prng = SystemRandom()

#         # cls.order = bn256.order
#         # cls._g1 = bn256.curve_G
#         # cls._g2 = bn256.twist_G
#         # cls.gf = GF(modulus = cls.order)
#         # cls.gf.is_signed = False
#         group1 = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
#         group2 = EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
#         group1.is_additive = False
#         group1.is_multiplicative = True
#         group2.is_additive = False
#         group2.is_multiplicative = True
#         cls.order = group1.order
#         cls._g1 = group1.generator
#         cls._g2 = group2.generator
#         cls.gf = GF(modulus = cls.order)

#     @classmethod
#     def tearDownClass(cls):
#         pass

#     def setUp(self):
#         # Randomness for commitment
#         self.gamma = self.gf(self.prng.randrange(1, self.order))


#     def test_open_linear_form_koe(self):
#         gf = self.gf

#         # Inputs
#         x = [gf(1), gf(2), gf(0), gf(0)]
#         n = len(x) 

#         # Trusted setup
#         pp = koe.trusted_setup(self._g1, self._g2, n, self.order)

#         # Create random linear form L
#         random_coeffs = list(gf(self.prng.randrange(self.order)) for i in range(n)) 
#         self.L = pivot.LinearForm(random_coeffs)

#         # Commit via [Gro10] restriction argument, the prove linear form opening
#         P, pi = koe.restriction_argument_prover(range(len(x)), x, self.gamma, pp)
#         proof, u = koe.opening_linear_form_prover(self.L, x, self.gamma, pp, P, pi)
#         verification = koe.opening_linear_form_verifier(self.L, pp, proof, u)
#         self.assertEqual(all(check == True for check in verification.values()), True)


#     def test_nullity_koe(self):
#         gf = self.gf

#         x = [gf(1), gf(2), gf(3)]
#         n = len(x)
#         lin_forms = [
#             pivot.LinearForm([6, 0, -2]), 
#             pivot.LinearForm([0, 3, -2]), 
#             pivot.LinearForm([2, 2, -2])]

#         # Trusted setup
#         pp = koe.trusted_setup(self._g1, self._g2, n, self.order)


#         """ Prover calculates proof for multiple linear forms 
#         using nullity protocol with KoE assumption.
#         """
#         S = range(len(x))
#         P, pi = koe.restriction_argument_prover(S, x, self.gamma, pp)
#         proof, L, u = koe.prove_nullity_koe(pp, lin_forms, x, self.gamma, gf, P, pi)

#         # Verifier checks proof
#         verification = koe.opening_linear_form_verifier(L, pp, proof, u)
#         self.assertEqual(all(check == True for check in verification.values()), True)


# if __name__ == '__main__':
#     unittest.main()