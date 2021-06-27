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

from mpyc.runtime import mpc
from mpyc.finfields import GF, PrimeFieldElement
from mpyc.sectypes import SecureFiniteField, SecureInteger

import verifiable_mpc.ac20_circuit_sat.circuit_builder as cb
import verifiable_mpc.ac20_circuit_sat.circuit_sat_cb as cs
import verifiable_mpc.ac20_circuit_sat.mpc_ac20_cb as mpc_cs
from sec_groups.fingroups import QuadraticResidue, EllipticCurve
import sec_groups.ellcurves as ell
from sec_groups.secgroups import SecureGroup, secure_repeat
from sec_groups.tools.find_primes import find_safe_primes


prng = SystemRandom()
pp = pprint.PrettyPrinter(indent=4)

PIVOT = (
    cs.PivotChoice.compressed
)  # pivot, compressed; TODO: koe is not supported, to be added (with BN256 curves)
GROUP = "QR"


async def main(pivot_choice, group_choice, n):
    # General setup
    if group_choice == "Elliptic":
        group = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
        group.is_additive = False
        group.is_multiplicative = True

        await mpc.start()
        sec_grp = SecureGroup(group, group.arithm)
        secint = mpc.SecInt(p=sec_grp.group.order)
        gf = secint.field
    elif group_choice == "QR":
        order, modulus = find_safe_primes(64)
        group = QuadraticResidue(modulus=modulus)

        await mpc.start()
        sec_grp = SecureGroup(group)
        secint = mpc.SecInt(p=sec_grp.group.order)
        gf = secint.field
    else:
        raise ValueError

    assert sec_grp.group.is_multiplicative
    print(f"Group: {GROUP}")

    circuit = cb.Circuit()
    b = cb.CircuitVar(secint(1), circuit, "b")
    c = cb.CircuitVar(secint(2), circuit, "c")
    d = c + c + c * c + c * c * 1 + 1 + b
    e = d * d + c + 10
    f = d * c + e
    f.label_output("f")
    g = f + 100
    g.label_output("g")
    h = g >= 10
    h.label_output("h")

    x = circuit.initial_inputs()
    # Check if resulting commitment vector is of appropriate length.
    check, padding, g_length = cs.check_input_length_power_of_2(x, circuit)
    # Add unused variables to pad the length of the commitment vector to power of 2 minus 1.
    dummies = [cb.CircuitVar(0, circuit, "unused_"+str(i)) for i in range(padding)]
    x = circuit.initial_inputs()
    print("Length of input vector including auxiliary inputs (witnesses for special gates): ", len(x))
    print("Length of commitment vector: ", g_length)

    print("Create generators.")
    generators = await mpc_cs.create_generators(group, secint, g_length)

    print("Start non-interactive circuit satisfiability proof with compressed pivot. ")
    proof = await mpc_cs.circuit_sat_prover(generators, circuit, x, gf, pivot_choice)

    print("Start verification.")
    verification = cs.circuit_sat_verifier(proof, generators, circuit, gf, pivot_choice)

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
    parser.add_argument("-n", type=int, help="roughly number of multiplications")
    parser.add_argument(
        "--elliptic",
        action="store_true",
        help="use elliptic curve groups (default QR groups)",
    )
    parser.set_defaults(n=100)
    args = parser.parse_args()
    if args.elliptic:
        GROUP = "Elliptic"
    verification = mpc.run(main(PIVOT, GROUP, args.n))


# TODO-list
# TODO-1: Reduce is_prime checks if redundant (for QR group)
# TODO-1: Add koe-assumption to verifiable MPC demo, using BN256 curves.
# TODO-2: Check if necessary to open z_prime in protocol_4_prover
