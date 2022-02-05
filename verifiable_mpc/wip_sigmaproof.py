# Work in progress
# Implement example(s) from exercise 5.3.2 from Cryptographic Protocols lecture 
# See page 53 of: https://www.win.tue.nl/~berry/CryptographicProtocols/LectureNotes.pdf 
# Goal: Combine shorter sigma-proofs for gadgets with circuit proofs

from random import SystemRandom
from mpyc.finfields import GF
from mpyc.fingroups import QuadraticResidues, EllipticCurve
import verifiable_mpc.ac20.pivot as pivot


# Setup
prng = SystemRandom()
group = QuadraticResidues(l=64)
g = group.generator
n = group.order
gf = GF(modulus=n)
r = gf(prng.randrange(1, n))
h = g^int(r)

# Prove x != 0. (Exercise 5.3.2.g., answer on p. 101)
# Relation {(B;x,y): B=g^x h^y, psi(x,y)}
# psi(x,y) = x!=0

# Prover
# B = g^x h^y, where x, y are private input to the prover
x, y = 1, 1
B = (g^x)*(h^y)
u, v = [gf(prng.randrange(1, n)) for i in range(2)]
a = (B^int(u))*(h^int(v))

# In interactive protocol, verifier samples c: c = gf(prng.randrange(1, n))
# Make non-interactive with Fiat-Shamir heuristic
input_list = [a, B]
c = gf(pivot.fiat_shamir_hash(input_list, gf.order))

# Prover
r = u + c/x
s = v - c*y/x

# Verifier
print((B^int(r))*(h^int(s)) == a*(g^int(c)))

# Secure groups
from mpyc.runtime import mpc
from mpyc.secgroups import repeat_public_base_public_output as secure_repeat
mpc.run(mpc.start())
secgrp = mpc.SecGrp(group)
sectype = mpc.SecFld(modulus=n)

# Setup
r = mpc._random(sectype)
h = mpc.run(secgrp.repeat_public(g, r))

# Prover
x, y = sectype(1), sectype(1)
B = mpc.run(secgrp.repeat_public([g, h], [x, y]))
u, v = [mpc._random(sectype) for i in range(2)]
a = mpc.run(secgrp.repeat_public([B, h], [u, v]))
input_list = [a, B]
c = gf(pivot.fiat_shamir_hash(input_list, gf.order))
r = mpc.run(mpc.gather(u + c/x))
s = mpc.run(mpc.gather(v - c*y/x))

# Verifier
print((B^int(r))*(h^int(s)) == a*(g^int(c)))

mpc.run(mpc.shutdown())
