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

import pprint
import os
import sys
from mpyc.fingroups import QuadraticResidues, EllipticCurve
from mpyc.runtime import mpc
import verifiable_mpc.ac20.circuit_builder as cb
import verifiable_mpc.ac20.circuit_sat_cb as cs
import verifiable_mpc.ac20.mpc_ac20_cb as mpc_cs


pp = pprint.PrettyPrinter(indent=4)

PIVOT = cs.PivotChoice.compressed  
GROUP = "QR"


async def main(pivot_choice, group_choice, n):
    await mpc.start()

    # General setup
    if pivot_choice == cs.PivotChoice.koe:
        group1 = EllipticCurve('BN256', 'projective')  # 'jacobian'
        group2 = EllipticCurve('BN256_twist', 'projective')  # 'jacobian'
        # NB: using projective here because we need oblivious implementation of groups
        group1.is_additive = False
        group1.is_multiplicative = True
        group2.is_additive = False
        group2.is_multiplicative = True
        group = [group1, group2]
        sec_grp = mpc.SecGrp(group1)
        sec_grp2 = mpc.SecGrp(group2)
        assert sec_grp.group.order == sec_grp2.group.order
    elif group_choice == "Elliptic":
        group = EllipticCurve('Ed25519', 'projective')
        group.is_additive = False
        group.is_multiplicative = True
        sec_grp = mpc.SecGrp(group)
    elif group_choice == "QR":
        group = QuadraticResidues(l=1024)
        sec_grp = mpc.SecGrp(group)
    else:
        raise ValueError
    
    print("Start AC20 with group: ", group)

    sectype = mpc.SecInt(l=16, p=sec_grp.group.order)
    print("sectype.bit_length:", sectype.bit_length)
    # sectype = mpc.SecFld(modulus=sec_grp.group.order)
    gf = sectype.field    

    circuit = cb.Circuit()
    # For integers, ensure enough (30 bit) headroom vs SecInt.p
    b = cb.CircuitVar(sectype(1), circuit, "b")
    c = cb.CircuitVar(sectype(2), circuit, "c")

    d = c + c + c * c + c * c * 1 + 1 + b
    e = d * d + c**n + 10
    f = d * c + e
    f.label_output("f")
    g = f != 100
    g.label_output("g")
    h = g >= 10  # Comparison requires sectype = mpc.Secint(..). Slow with KoE pivot.
    h.label_output("h")

    x = circuit.initial_inputs()
    # Check if resulting commitment vector is of appropriate length.
    check, padding, g_length = cs.check_input_length_power_of_2(x, circuit)
    # Add unused variables to pad the length of the commitment vector to power of 2 minus 1.
    unused = [cb.CircuitVar(0, circuit, "unused_"+str(i)) for i in range(padding)]
    x = circuit.initial_inputs()

    print("Length of input vector including auxiliary inputs (witnesses for special gates): ", len(x))
    print("Length of commitment vector: ", g_length)

    print("Create generators.")
    if pivot_choice in [cs.PivotChoice.pivot, cs.PivotChoice.compressed]:
        generators = await mpc_cs.create_generators(group, sectype, g_length)
    elif pivot_choice in [cs.PivotChoice.koe]:
        generators = await mpc_cs.koe_trusted_setup(group, sectype, g_length, progress_bar=True)
    else:
        raise NotImplementedError

    print("Start non-interactive circuit satisfiability proof. ")
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
    parser.add_argument("--elliptic", action="store_true", 
                        help="use elliptic curve groups (default QR groups)")
    parser.add_argument('--basic', action='store_true',
                        help='use basic pivot (not the compressed pivot)')
    parser.add_argument('--koe', action='store_true',
                        help='use pivot based on Knowledge-of-Exponent assumption and BN256 curves')
    parser.set_defaults(n=3)
    args = parser.parse_args()
    if args.elliptic:
        GROUP = "Elliptic"
    elif args.basic:
        PIVOT = cs.PivotChoice.pivot
    elif args.koe:
        PIVOT = cs.PivotChoice.koe

    verification = mpc.run(main(PIVOT, GROUP, args.n))


# TODO-list
# TODO-2: improve syntax for passing two groups to create_generators
# TODO-2: Check if necessary to open z_prime in protocol_4_prover
