import sys
import os
import unittest
from random import SystemRandom
import verifiable_mpc.ac20.pivot as pivot
from mpyc.finfields import GF
from mpyc.fingroups import QuadraticResidues


class Pivot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prng = SystemRandom()
        cls.group = QuadraticResidues(l=64)
        cls.order = cls.group.order
        cls.h = cls.group.generator
        cls.gf = GF(modulus=cls.order)

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        # Define input
        gf = self.gf
        self.x = [gf(1), gf(2), gf(0), gf(0)]

        # Create random generators
        self.n = len(self.x)
        random_exponents = []
        while len(set(random_exponents)) != self.n:
            random_exponents = list(
                self.prng.randrange(1, self.order) for i in range(self.n)
            )
        self.g = [self.h ** i for i in random_exponents]
        self.assertEqual(len({int(g_i) for g_i in self.g}), self.n)

    def test_pivot_interactive(self):
        # Create random element for commitments
        gamma = self.prng.randrange(self.order)
        # Public parameters
        P = pivot.vector_commitment(self.x, gamma, self.g, self.h)
        # Create random linear form L for testing
        random_coeffs = list(self.prng.randrange(self.order) for i in range(self.n))
        L = pivot.LinearForm(random_coeffs)
        y = L(self.x)

        # Interactive proof
        # Step 1: Prover announces with t, A
        r = list(self.prng.randrange(self.order) for i in range(self.n))
        rho = self.prng.randrange(self.order)
        t = L(r)
        A = pivot.vector_commitment(r, rho, self.g, self.h)

        # Step 2: Verifier responds with challenge
        c = self.prng.randrange(self.order)

        # Step 3: Prover responds
        z = [c * x_i + r[i] for i, x_i in enumerate(self.x)]
        phi = c * gamma + rho

        # Step 4: Verifier checks
        self.assertEqual(pivot.vector_commitment(z, phi, self.g, self.h), A * (P ** c))
        self.assertEqual(L(z), c * y + t)

    def test_pivot_noninteractive(self):
        # Create random element for commitments
        gamma = self.prng.randrange(self.order)
        # Public parameters
        P = pivot.vector_commitment(self.x, gamma, self.g, self.h)
        # Create random linear form L for testing
        random_coeffs = list(self.prng.randrange(self.order) for i in range(self.n))
        L = pivot.LinearForm(random_coeffs)
        y = L(self.x)

        z, phi, c = pivot.prove_linear_form_eval(
            self.g, self.h, P, L, y, self.x, gamma, self.gf
        )
        verification = pivot.verify_linear_form_proof(
            self.g, self.h, P, L, y, z, phi, c
        )
        self.assertEqual(verification, True)

    def test_linear_form(self):
        # Test LinearForm and AffineForm classes
        lf = pivot.LinearForm([0, 1, 2])
        self.assertEqual(
            (lf + lf + 2 * lf + lf.eval([1, 1, 1]) - lf).eval([1, 2, 3]), 27
        )
        self.assertEqual(lf([1, 2, 3]), 8)


if __name__ == "__main__":
    unittest.main()
