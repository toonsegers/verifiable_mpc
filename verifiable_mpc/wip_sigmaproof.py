# Work in progress
# Implement example(s) from exercise 5.3.2 from Cryptographic Protocols lecture 
# See page 53 of: https://www.win.tue.nl/~berry/CryptographicProtocols/LectureNotes.pdf 
# Goal: Combine shorter sigma-proofs for gadgets with circuit proofs

from enum import Enum
from random import SystemRandom
from mpyc.finfields import GF
from mpyc.fingroups import QuadraticResidues, EllipticCurve
from mpyc.runtime import mpc
from mpyc.secgroups import repeat_public_base_public_output as secure_repeat
import verifiable_mpc.ac20.pivot as pivot
from mpyc.sectypes import SecureFiniteField, SecureInteger


class SigmaProof(Enum):
    """Predicate of proof. """

    not_zero = 1


prng = SystemRandom()


async def sigma_prove_not_zero(x, group):
    """Generate sigma proof that x != 0, using the Discrete Log assumption for given group.

    Prove x != 0. (Exercise 5.3.2.g of Cryptographic Protocols, answer on p. 101)
    Relation {(B;x,y): B=g^x h^y, psi(x,y)}
    psi(x,y) = x!=0
    """
    # TODO: Consider separating generators and commitment from return value
    gf = GF(modulus=group.order)
    g = group.generator
    y = 1  # TODO: Redundant?

    if isinstance(x, (SecureFiniteField, SecureInteger)):
        sectype = type(x)
        secgrp = mpc.SecGrp(group)
        r = mpc._random(sectype)
        h = await secgrp.repeat_public(g, r)

        # Prover
        B = await secgrp.repeat_public([g, h], [x, y])
        u, v = [mpc._random(sectype) for i in range(2)]
        a = await secgrp.repeat_public([B, h], [u, v])
        input_list = [a, B]
        c = gf(pivot.fiat_shamir_hash(input_list, gf.order))
        r = await mpc.gather(u + c/x)
        s = await mpc.gather(v - c*y/x)
    elif isinstance(x, int):
        r = gf(prng.randrange(1, group.order))
        h = g^int(r)
        # Prover
        # B = g^x h^y, where x, y are private input to the prover
        B = (g^x)*(h^y)
        u, v = [gf(prng.randrange(1, group.order)) for i in range(2)]
        a = (B^int(u))*(h^int(v))
        # In interactive protocol, verifier samples c: c = gf(prng.randrange(1, n))
        # Make non-interactive with Fiat-Shamir heuristic
        input_list = [a, B]
        c = gf(pivot.fiat_shamir_hash(input_list, gf.order))
        r = u + c/x
        s = v - c*y/x
    else:
        raise TypeError

    proof = {
        "predicate": SigmaProof.not_zero,
        "generators": (g, h),
        "commitment": B,
        "proof": (a, r, s),
        }
    return proof


def sigma_verify_not_zero(proof):
    assert proof["predicate"] == SigmaProof.not_zero
    g, h = proof["generators"]
    B = proof["commitment"]
    a, r, s = proof["proof"]
    input_list = [a, B]
    gf = GF(modulus=type(g).order)
    c = gf(pivot.fiat_shamir_hash(input_list, gf.order))
    return (B^int(r))*(h^int(s)) == a*(g^int(c))


async def example_not_zero():
    await mpc.start()
    group = QuadraticResidues(l=64)
    secgrp = mpc.SecGrp(group)
    sectype = mpc.SecFld(modulus=group.order)
    # Prover
    # x = sectype(1)
    x = 1
    proof = await sigma_prove_not_zero(x, group)
    # Verifier
    print(sigma_verify_not_zero(proof))
    await mpc.shutdown()


if __name__ == "__main__":
    mpc.run(example_not_zero())


