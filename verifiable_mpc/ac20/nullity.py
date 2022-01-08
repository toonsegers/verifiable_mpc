""" Implementation of https://eprint.iacr.org/2020/152
``Compressed Σ-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

Protocols:
* Protocol Pi_Nullity, page 17-18, to proof null-evaluations with "polynomial amortization trick" using pivot

"""

import os
import sys
from random import SystemRandom
import pprint
import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.compressed_pivot as compressed_pivot

prng = SystemRandom()
pp = pprint.PrettyPrinter(indent=4)


def prove_nullity_compressed(generators, P, lin_forms, x, gamma, gf):

    input_list = [P, lin_forms]
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    y = L(x)
    proof = compressed_pivot.protocol_5_prover(generators, P, L, y, x, gamma, gf)
    return proof, L, y, rho


def verify_nullity_compressed(generators, P, L, lin_forms, rho, y, proof, gf):
    L_check = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    if not L_check == L:
        print(
            "Linear form L does not correspond to reconstructed linear form with rho."
        )
        verification = False
        return verification
    verification = compressed_pivot.protocol_5_verifier(generators, P, L, y, proof, gf)
    return verification
