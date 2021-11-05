"""Implementation of Trinocchio [SVV15] in Python.

Paper by Schoenmakers, Veeningen and De Vreede:
https://eprint.iacr.org/2015/480

Adaptations to original Trinocchio protocol:
* TBD

Credits:
* Meilof Veeningen's PySNARK: https://github.com/meilof/pysnark
"""

import os, sys

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.runtime import mpc
from mpyc.fingroups import EllipticCurve
from mpyc.fingroups import FiniteGroupElement


bn_curve = EllipticCurve('BN256', 'jacobian')
g1 = bn_curve.generator
bn_twist = EllipticCurve('BN256_twist', 'jacobian')
g2 = bn_twist.generator
modulus = bn_curve.order
point_add = FiniteGroupElement.__matmul__


# TODO list
# TODO-1: Implement ZK variant
# TODO-2: Remove progress bar (Bar)
# TODO-1: Apply Mark's trick consistently
# TODO-1: Implement multi-client with basic Trinocchio
# TODO-1: Implement full Trinocchio scheme (Alg. 4, p. 25)
# TODO-1: Enable __neg__ and __sub__ in bn256 classes, and then don't set gf.is_signed=False
# TODO-3: Loosely couple/decouple type detection from "g1"/"g2" in key name in: pt_type = bn256.CurvePoint if "g1" in key else bn256.CurveTwist
