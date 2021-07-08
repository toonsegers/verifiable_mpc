"""Demo ``Compressed Sigma-Protocol Theory'' by Attema and Cramer. 

Full title: ``Compressed Sigma-Protocol Theory and Practical 
Application to Plug & Play Secure Algorithmics'' 
by Thomas Attema and Ronald Cramer.
Link: https://eprint.iacr.org/2020/152

This demo proves circuit satisfiability using the pivots
explained in the paper (regular, compressed and KoE-based pivot).

Credits:
* r1cs and qap tools by Vitalik Buterin:
https://github.com/ethereum/research/tree/master/zksnark (MIT license)
(Adapted for this purpose).
"""

from random import SystemRandom
import pprint
import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.finfields import GF

import verifiable_mpc.ac20.circuit_sat_cb as cs
import verifiable_mpc.ac20.circuit_builder as cb
from sec_groups.fingroups import QuadraticResidue, EllipticCurve
import sec_groups.ellcurves as ell
from sec_groups.tools.find_primes import find_safe_primes


pp = pprint.PrettyPrinter(indent=4)


GROUP = "QR"  # "Elliptic" or "QR"


def main(pivot_choice):
    print("Pivot selected: ", pivot_choice)

    if pivot_choice == cs.PivotChoice.koe:
        # TODO: improve syntax for passing two groups to create_generators
        group1 = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        group2 = EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        group1.is_additive = False
        group1.is_multiplicative = True
        group2.is_additive = False
        group2.is_multiplicative = True
        group = [group1, group2]
        order = group1.order
        gf = GF(modulus=order)
    elif GROUP == "Elliptic":
        group = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
        group.is_additive = False
        group.is_multiplicative = True
        gf = GF(modulus=group.order)
    elif GROUP == "QR":
        order, modulus = find_safe_primes(64)
        group = QuadraticResidue(modulus=modulus)
        gf = GF(modulus=group.order)

    # Inputs: Python code to evaluate and x-vector.
    circuit = cb.Circuit()
    b = cb.CircuitVar(gf(1), circuit, "b")
    c = cb.CircuitVar(gf(2), circuit, "c")

    d = c + c + c * c + c * c * 1 + 1 + b 
    e = d*d + c + 10
    f = d*c + e
    f.label_output("f")
    print("f=", f)
    g = f + 100
    g.label_output("g")
    print("g=", g)
    h = g != 10
    h.label_output("h")

    x = circuit.initial_inputs()
    # Check if resulting commitment vector is of appropriate length.
    check, padding, g_length = cs.check_input_length_power_of_2(x, circuit)
    # Add unused variables to pad the length of the commitment vector to power of 2 minus 1.
    dummies = [cb.CircuitVar(0, circuit, "unused_"+str(i)) for i in range(padding)]
    x = circuit.initial_inputs()
    print("Length of input vector including auxiliary inputs (witnesses for special gates): ", len(x))
    print("Length of commitment vector: ", g_length)

    generators = cs.create_generators(g_length, pivot_choice, group, progress_bar=True)
    print("Generators created/trusted setup done.")

    print("Start non-interactive circuit satisfiability proof with compressed pivot. ")
    proof = cs.circuit_sat_prover(generators, circuit, x, gf, pivot_choice)
    # print("Proof:")
    # pp.pprint(proof)
    print("Start verification.")
    verification = cs.circuit_sat_verifier(proof, generators, circuit, gf, pivot_choice)
    print("Verification checks: ")
    pp.pprint(verification)

    return verification



if __name__ == "__main__":
    pivot_options = {
        1: cs.PivotChoice.pivot,
        2: cs.PivotChoice.compressed,
        3: cs.PivotChoice.koe,
    }
    choice = int(input(f"Select {pivot_options}, default [2]: ") or "2")
    pivot_choice = pivot_options[choice]

    verification = main(pivot_choice)



# TODO:
# TS-1: Make interface of create_generators(,,group,) more intuitive when working with KOE
# TS-1: Remove L from proof; non-interactive nullity proof does not have to pass it to verifier.
