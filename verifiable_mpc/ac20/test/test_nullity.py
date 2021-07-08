import sys, os

project_root = sys.path.append(os.path.abspath("../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from random import SystemRandom
import unittest

import verifiable_mpc.ac20.circuit_sat_r1cs as cs
import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.nullity as nullity
import verifiable_mpc.ac20.compressed_pivot as compressed_pivot
from mpyc.finfields import GF, FiniteFieldElement
from sec_groups.tools.find_primes import find_safe_primes
from sec_groups.fingroups import QuadraticResidue


class Nullity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prng = SystemRandom()

        cls.order, cls.modulus = find_safe_primes(64)
        cls.group = QuadraticResidue(modulus=cls.modulus)
        cls.gf = GF(modulus=cls.group.order)

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        # Define input
        gf = self.gf
        self.x = [gf(1), gf(2), gf(3)]
        self.lin_forms = [
            pivot.LinearForm([6, 0, -2]),
            pivot.LinearForm([0, 3, -2]),
            pivot.LinearForm([2, 2, -2]),
        ]

        # Create random generators
        self.n = len(self.x)
        self.generators = cs.create_generators(
            self.n, cs.PivotChoice.compressed, self.group, progress_bar=False
        )
        self.g = self.generators["g"]
        self.h = self.generators["h"]
        self.k = self.generators["k"]

        self.s = len(self.lin_forms)
        # Create random element for commitments
        self.gamma = self.prng.randrange(self.order)
        # Verifier samples random rho
        self.rho = self.prng.randrange(self.order)

        # Prover to open L(x) = Sum L_i(x)*rho**(i-1) for i = 1 .. s
        self.P = pivot.vector_commitment(
            self.x, self.gamma, self.generators["g"], self.generators["h"]
        )
        self.L = sum((self.lin_forms[i]) * (self.rho ** i) for i in range(self.s))
        self.y = self.L(self.x)

    def test_nullity_compressed_noninteractive(self):
        # Non-interactive version of nullity check

        proof, L, y, rho = nullity.prove_nullity_compressed(
            self.generators, self.P, self.lin_forms, self.x, self.gamma, self.gf
        )
        verification = nullity.verify_nullity_compressed(
            self.generators, self.P, L, self.lin_forms, rho, y, proof, self.gf
        )
        self.assertEqual(verification, True)


if __name__ == "__main__":
    unittest.main()
