"""Demo prover side of Circuit Sat ZK from [AC20] in MPC.

Based on: "Compressed Sigma-Protocol Theory and Practical
Application to Plug & Play Secure Algorithmics"
by Thomas Attema and Ronald Cramer.
Link: https://eprint.iacr.org/2020/152

Uses MPyC framework and 'Secure Groups' abstraction.

Credits:
* Verifiable MPC introduced by Berry Schoenmakers, Meilof
  Veeningen and Niels de Vreede: https://eprint.iacr.org/2015/480
* "Compressed Sigma-Protocol Theory" by Thomas Attema and Ronald Cramer.
"""

from random import SystemRandom
import pprint
import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import verifiable_mpc.ac20_circuit_sat.circuit_sat_r1cs as cs
import verifiable_mpc.ac20_circuit_sat.mpc_ac20 as mpc_cs
from mpyc.runtime import mpc
from mpyc.finfields import GF, PrimeFieldElement
from mpyc.sectypes import SecureFiniteField, SecureInteger
from sec_groups.tools.find_primes import find_safe_primes
from sec_groups.fingroups import QuadraticResidue, EllipticCurve
import sec_groups.ellcurves as ell
from sec_groups.secgroups import SecureGroup, secure_repeat


prng = SystemRandom()
pp = pprint.PrettyPrinter(indent=4)

# PIVOT = cs.PivotChoice.compressed   # pivot, compressed; TODO: koe is not supported, to be added (with BN256 curves)
# GROUP = "QR"
PIVOT = cs.PivotChoice.koe   # todo
GROUP = "KoE"


async def main(pivot_choice, group_choice, n):
    await mpc.start()

    # General setup
    if group_choice == "Elliptic":
        group = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
        group.is_additive = False
        group.is_multiplicative = True

        sec_grp = SecureGroup(group, group.arithm)
        secfld = mpc.SecFld(modulus=sec_grp.group.order)
        gf = secfld.field
    elif group_choice == "QR":
        order, modulus = find_safe_primes(2048)
        group = QuadraticResidue(modulus=modulus)

        sec_grp = SecureGroup(group)
        secfld = mpc.SecFld(modulus=sec_grp.group.order)
        gf = secfld.field
    elif group_choice == "KoE":
        # TODO: improve syntax for passing two groups to create_generators
        group1 = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        group2 = EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        group1.is_additive = False
        group1.is_multiplicative = True
        group2.is_additive = False
        group2.is_multiplicative = True
        group = [group1, group2]
        order = group1.order

        sec_grp1 = SecureGroup(group1)
        sec_grp2 = SecureGroup(group2)
        assert sec_grp1.group.order == sec_grp2.group.order
        secfld = mpc.SecFld(modulus=sec_grp1.group.order)
        gf = secfld.field
    else:
        raise ValueError

    # assert sec_grp.group.is_multiplicative

    print(f"Group: {GROUP}")

    code = f"""
def qeval(x1, x2):
    z = 2*x1 + 3*x2
    y = x1**{n}
    return y
"""
    x = [secfld(1), secfld(2)]
    x, code, g_length = cs.input_length_power_of_2(x, code)
    print(
        "Code to prove circuit suitasfiability for (after padding to ensure length(z) + 1 is power of 2):"
    )
    print(code)
    print("Length of commitment vector=", g_length)
    print("Create generators.")
    if pivot_choice in [cs.PivotChoice.pivot, cs.PivotChoice.compressed]:
        generators = await mpc_cs.create_generators(group, secfld, g_length)
    elif pivot_choice in [cs.PivotChoice.koe]:
        generators = await mpc_cs.koe_trusted_setup(group, secfld, g_length, progress_bar=True)
    else:
        raise NotImplementedError

    print("Start non-interactive circuit satisfiability proof with compressed pivot. ")
    proof = await mpc_cs.circuit_sat_prover(generators, code, x, gf, pivot_choice)

    print("Start verification.")
    verification = cs.circuit_sat_verifier(proof, generators, code, gf, pivot_choice)

    if all(verification.values()):
        print("All checks passed.")
    else:
        print("Not all checks passed.")
    pp.pprint(verification)

    await mpc.shutdown()
    return verification


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int,
                        help='roughly number of multiplications')
    parser.add_argument('--elliptic', action='store_true',
                        help='use elliptic curve groups (default QR groups)')
    parser.set_defaults(n=10)
    args = parser.parse_args()
    if args.elliptic:
        GROUP = "Elliptic"
    verification = mpc.run(main(PIVOT, GROUP, args.n))


# TODO-list
# TODO-1: Reduce is_prime checks if redundant (for QR group)
# TODO-1: Add koe-assumption to verifiable MPC demo, using BN256 curves.
# TODO-2: Check if necessary to open z_prime in protocol_4_prover
