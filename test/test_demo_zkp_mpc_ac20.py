import sys, os
import unittest
import verifiable_mpc.ac20.circuit_sat_cb as cs
from demos import demo_zkp_mpc_ac20
from mpyc.runtime import mpc


class CircuitSat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_cs_pivot_elliptic(self):
        verification = mpc.run(demo_zkp_mpc_ac20.main(cs.PivotChoice.pivot, "Elliptic", 2))
        self.assertEqual(all(check == True for check in verification.values()), True)

    def test_cs_compressed_elliptic(self):
        verification = mpc.run(demo_zkp_mpc_ac20.main(cs.PivotChoice.compressed, "Elliptic", 2))
        self.assertEqual(all(check == True for check in verification.values()), True)

    def test_cs_knowledgeofexponent_elliptic(self):
        verification = mpc.run(demo_zkp_mpc_ac20.main(cs.PivotChoice.koe, "Elliptic", 2))
        self.assertEqual(verification["y1*y2=y3"], True)
        self.assertEqual(verification["L_wellformed_from_Cfgh_forms"], True)
        self.assertEqual(verification["pivot_verification"]["restriction_arg_check"], True)
        self.assertEqual(verification["pivot_verification"]["PRQ_check"], True)

    def test_cs_pivot_quadratic_residues(self):
        verification = mpc.run(demo_zkp_mpc_ac20.main(cs.PivotChoice.pivot, "QR", 10))
        self.assertEqual(all(check == True for check in verification.values()), True)

    def test_cs_compressed_quadratic_residues(self):
        verification = mpc.run(demo_zkp_mpc_ac20.main(cs.PivotChoice.compressed, "QR", 10))
        self.assertEqual(all(check == True for check in verification.values()), True)


if __name__ == "__main__":
    unittest.main()
