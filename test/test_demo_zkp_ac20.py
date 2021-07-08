import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import unittest

import verifiable_mpc.ac20.circuit_sat_cb as cs
from demos import demo_zkp_ac20


class CircuitSat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_cs_pivot(self):
        verification = demo_zkp_ac20.main(cs.PivotChoice.pivot)
        self.assertEqual(verification["y1*y2=y3"], True)
        self.assertEqual(verification["L_wellformed_from_Cfgh_forms"], True)
        self.assertEqual(verification["pivot_verification"], True)

    def test_cs_compressed(self):
        verification = demo_zkp_ac20.main(cs.PivotChoice.compressed)
        self.assertEqual(verification["y1*y2=y3"], True)
        self.assertEqual(verification["L_wellformed_from_Cfgh_forms"], True)
        self.assertEqual(verification["pivot_verification"], True)

    def test_cs_koe(self):
        verification = demo_zkp_ac20.main(cs.PivotChoice.koe)
        self.assertEqual(verification["y1*y2=y3"], True)
        self.assertEqual(verification["L_wellformed_from_Cfgh_forms"], True)
        self.assertEqual(
            verification["pivot_verification"]["restriction_arg_check"], True
        )
        self.assertEqual(verification["pivot_verification"]["PRQ_check"], True)


if __name__ == "__main__":
    unittest.main()
