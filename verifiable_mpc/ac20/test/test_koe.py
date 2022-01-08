import sys
import os
import unittest
from random import SystemRandom
import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.knowledge_of_exponent as koe
from mpyc.finfields import GF
from mpyc.fingroups import QuadraticResidues, EllipticCurve


class KoE(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.prng = SystemRandom()
        group1 = EllipticCurve('BN256', 'projective')  # 'jacobian'
        group2 = EllipticCurve('BN256_twist', 'projective')  # 'jacobian'
        group1.is_additive = False
        group1.is_multiplicative = True
        group2.is_additive = False
        group2.is_multiplicative = True
        cls.order = group1.order
        cls._g1 = group1.generator
        cls._g2 = group2.generator
        cls.gf = GF(modulus = cls.order)

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        # Randomness for commitment
        self.gamma = self.gf(self.prng.randrange(1, self.order))


    def test_open_linear_form_koe(self):
        gf = self.gf

        # Inputs
        x = [gf(1), gf(2), gf(0), gf(0)]
        n = len(x)

        # Trusted setup
        pp = koe.trusted_setup(self._g1, self._g2, n, self.order)

        # Create random linear form L
        random_coeffs = list(gf(self.prng.randrange(self.order)) for i in range(n))
        self.L = pivot.LinearForm(random_coeffs)

        # Commit via [Gro10] restriction argument, the prove linear form opening
        P, pi = koe.restriction_argument_prover(range(len(x)), x, self.gamma, pp)
        proof, u = koe.opening_linear_form_prover(self.L, x, self.gamma, pp, P, pi)
        verification = koe.opening_linear_form_verifier(self.L, pp, proof, u)
        self.assertEqual(all(check == True for check in verification.values()), True)


    def test_nullity_koe(self):
        gf = self.gf

        x = [gf(1), gf(2), gf(3)]
        n = len(x)
        lin_forms = [
            pivot.LinearForm([6, 0, -2]),
            pivot.LinearForm([0, 3, -2]),
            pivot.LinearForm([2, 2, -2])]

        # Trusted setup
        pp = koe.trusted_setup(self._g1, self._g2, n, self.order)


        """ Prover calculates proof for multiple linear forms
        using nullity protocol with KoE assumption.
        """
        S = range(len(x))
        P, pi = koe.restriction_argument_prover(S, x, self.gamma, pp)
        proof, L, u = koe.prove_nullity_koe(pp, lin_forms, x, self.gamma, gf, P, pi)

        # Verifier checks proof
        verification = koe.opening_linear_form_verifier(L, pp, proof, u)
        self.assertEqual(all(check == True for check in verification.values()), True)


if __name__ == '__main__':
    unittest.main()