""" Adaptation of work by Jack Lloyd:
https://github.com/randombit/pairings.py/

Original code under BSD-2-Clause license:
https://github.com/randombit/pairings.py/blob/master/LICENSE.md

----------------------------------------------------------------------
Copyright and license information of the original work.
----------------------------------------------------------------------

Pairings over a 256-bit BN curve
(C) 2017 Jack Lloyd <jack@randombit.net>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

----------------------------------------------------------------------
Note: Code below is adapted from its original to directly use extension
fields as defined in MPyC (mpyc.finfields.ExtensionFieldElement).
----------------------------------------------------------------------

"""
import copy

from mpyc.fingroups import EllipticCurve
BN256 = EllipticCurve('BN256')
BN256_TWIST = EllipticCurve('BN256_twist')

GFp_1 = BN256.field
GFp_2 = BN256_TWIST.field

v = 1868033
u = pow(v, 3)
p = 36 * u ** 4 + 36 * u ** 3 + 24 * u ** 2 + 6 * u + 1

i = GFp_2([0, 1, 0])
xi = i + GFp_2([3, 0, 0])

xi1 = [
    xi ** (1 * (p - 1) // 6),
    xi ** (2 * (p - 1) // 6),
    xi ** (3 * (p - 1) // 6),
    xi ** (4 * (p - 1) // 6),
    xi ** (5 * (p - 1) // 6),
]


gfp_2_zero = GFp_2(0)
gfp_2_one = GFp_2(1)


def conjugate(a):
    assert isinstance(a, BN256_TWIST.field)
    a_const = a.value.value[0]
    a_i = a.value.value[1]
    return type(a)([a_const, p - a_i])


xi2 = [(x * conjugate(x)) for x in xi1]


def to_naf(x):
    z = []
    while x > 0:
        if x % 2 == 0:
            z.append(0)
        else:
            zi = 2 - (x % 4)
            x -= zi
            z.append(zi)
        x = x // 2
    return z


# 6u+2 in NAF
naf_6up2 = list(reversed(to_naf(6 * u + 2)))[1:]


def bits_of(k):
    return [int(c) for c in "{0:b}".format(k)]


# cubic extension of GFp_2
class GFp_6(object):
    def __init__(self, x, y, z):
        assert type(x) == GFp_2 and type(y) == GFp_2 and type(z) == GFp_2
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __repr__(self):
        return "(%s,%s,%s)" % (self.x, self.y, self.z)

    def is_zero(self):
        return int(self) == 0

    def is_one(self):
        return int(self) == 1

    def negative_of(self):
        return GFp_6(-self.x, -self.y, -self.z)

    def add(a, b):
        return GFp_6(a.x + b.x, a.y + b.y, a.z + b.z)

    def sub(a, b):
        return GFp_6(a.x - b.x, a.y - b.y, a.z - b.z)

    def double(self):
        return GFp_6(2 * self.x, 2 * self.y, 2 * self.z)

    def mul(a, b):
        # Algorithm 13 from http://eprint.iacr.org/2010/354.pdf
        # plus some short-circuits

        # if a.x.is_zero():
        # if a.y.is_zero():
        if int(a.x) == 0:
            if int(a.y) == 0:
                return b.mul_scalar(a.z)

            t0 = b.z * a.z
            t1 = b.y * a.y

            tz = (b.x + b.y) * (a.y)
            tz -= t1
            tz = tz * xi
            tz += t0

            ty = (b.y + b.z) * (a.y + a.z)
            ty -= t0
            ty -= t1

            tx = (b.x) * (a.z)
            tx += t1

            return GFp_6(tx, ty, tz)

        if int(b.x) == 0:
            if int(b.y) == 0:
                return a.mul_scalar(b.z)

            t0 = a.z * b.z
            t1 = a.y * b.y

            tz = (a.x + a.y) * (b.y)
            tz -= t1
            tz = tz * xi
            tz += t0

            ty = (a.y + a.z) * (b.y + b.z)
            ty -= t0
            ty -= t1

            tx = (a.x) * (b.z)
            tx += t1

            return GFp_6(tx, ty, tz)

        t0 = a.z * b.z
        t1 = a.y * b.y
        t2 = a.x * b.x

        tz = (a.x + a.y) * (b.x + b.y)
        tz -= t1
        tz -= t2
        tz = tz * xi
        tz += t0

        ty = (a.y + a.z) * (b.y + b.z)
        ty -= t0
        ty -= t1
        ty += t2 * xi

        tx = (a.x + a.z) * (b.x + b.z)
        tx -= t0
        tx += t1
        tx -= t2

        return GFp_6(tx, ty, tz)

    def __mul__(a, b):
        return a.mul(b)

    def __add__(a, b):
        return a.add(b)

    def __sub__(a, b):
        return a.sub(b)

    def mul_scalar(self, k):
        assert type(k) == GFp_2
        return GFp_6(self.x * k, self.y * k, self.z * k)

    def mul_tau(a):
        tx = a.y
        ty = a.z
        tz = a.x * xi
        return GFp_6(tx, ty, tz)

    def square(a):
        # Algorithm 16 from http://eprint.iacr.org/2010/354.pdf
        ay2 = a.y * 2
        c4 = a.z * ay2
        c5 = a.x ** 2
        c1 = c5 * xi + c4
        c2 = c4 - c5
        c3 = a.z ** 2
        c4 = a.x + a.z - a.y
        c5 = ay2 * a.x
        c4 = c4 ** 2
        c0 = c5 * xi + c3
        c2 = c2 + c4 + c5 - c3
        n = GFp_6(c2, c1, c0)
        return n

    def inverse(a):
        # Algorithm 17
        XX = a.x ** 2
        YY = a.y ** 2
        ZZ = a.z ** 2

        XY = a.x * a.y
        XZ = a.x * a.z
        YZ = a.y * a.z

        A = ZZ - XY * xi
        B = XX * xi - YZ
        # There is an error in the paper for this line
        C = YY - XZ

        F = (C * a.y) * xi
        F += A * a.z
        F += (B * a.x) * xi

        F = F.reciprocal()

        c_x = C * F
        c_y = B * F
        c_z = A * F
        return GFp_6(c_x, c_y, c_z)


gfp_6_zero = GFp_6(gfp_2_zero, gfp_2_zero, gfp_2_zero)
gfp_6_one = GFp_6(gfp_2_zero, gfp_2_zero, gfp_2_one)


class GFp_12(object):
    def __init__(self, x, y=None):
        assert type(x) == GFp_6
        assert type(y) == GFp_6
        self.x = x
        self.y = y

    def __mul__(self, other):
        return self.mul(other)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return "(%s,%s)" % (self.x, self.y)

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self):
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self):
        return GFp_12(self.x.negative_of(), self.y)

    def negative_of(self):
        return GFp_12(self.x.negative_of(), self.y.negative_of())

    def frobenius(self):
        e1_x = conjugate(self.x.x) * xi1[4]
        e1_y = conjugate(self.x.y) * xi1[2]
        e1_z = conjugate(self.x.z) * xi1[0]

        e2_x = conjugate(self.y.x) * xi1[3]
        e2_y = conjugate(self.y.y) * xi1[1]
        e2_z = conjugate(self.y.z)

        return GFp_12(GFp_6(e1_x, e1_y, e1_z), GFp_6(e2_x, e2_y, e2_z))

    def frobenius_p2(self):
        e1_x = self.x.x * xi2[4]
        e1_y = self.x.y * xi2[2]
        e1_z = self.x.z * xi2[0]

        e2_x = self.y.x * xi2[3]
        e2_y = self.y.y * xi2[1]
        e2_z = self.y.z

        return GFp_12(GFp_6(e1_x, e1_y, e1_z), GFp_6(e2_x, e2_y, e2_z))

    def sub(a, b):
        return GFp_12(a.x - b.x, a.y - b.y)

    def mul(a, b):
        # TODO Karatsuba (algo 20)
        AXBX = a.x * b.x
        AXBY = a.x * b.y
        AYBX = a.y * b.x
        AYBY = a.y * b.y
        return GFp_12(AXBY + AYBX, AYBY + AXBX.mul_tau())

    def mul_scalar(self, k):
        assert type(k) == GFp_6
        return GFp_12(self.x.mul(k), self.y.mul(k))

    def exp(self, k):
        assert isinstance(k, int)

        R = [GFp_12(gfp_6_zero, gfp_6_one), self]

        for kb in bits_of(k):
            R[kb ^ 1] = R[kb].mul(R[kb ^ 1])
            R[kb] = R[kb].square()

        return R[0]

    def square(a):
        v0 = a.x * a.y
        t = a.x.mul_tau()
        t += a.y
        ty = a.x + a.y
        ty *= t
        ty -= v0
        t = v0.mul_tau()
        ty -= t

        c_x = v0.double()
        c_y = ty

        return GFp_12(c_x, c_y)

    def inverse(a):
        e = GFp_12(a.x.negative_of(), a.y)

        t1 = a.x.square()
        t2 = a.y.square()
        t1 = t1.mul_tau()
        t1 = t2 - t1
        t2 = t1.inverse()

        e = e.mul_scalar(t2)
        return e


def line_func_add(r, p, q, r2):
    """ Requires twisted BN256 points r, p, BN256 point q and GFp_2 element r2."""

    assert type(r2) == GFp_2

    r_t = r.z ** 2
    B = p.x * r_t
    D = p.y + r.z
    D = D ** 2
    D -= r2
    D -= r_t
    D *= r_t

    H = B - r.x
    I = H ** 2

    E = I * 4

    J = H * E
    L1 = D - r.y
    L1 -= r.y

    V = r.x * E

    r_x = L1 ** 2
    r_x -= J
    r_x -= V * 2

    r_z = r.z + H
    r_z = r_z ** 2
    r_z -= r_t
    r_z -= I

    t = V - r_x
    t *= L1
    t2 = r.y * J
    t2 = t2 * 2
    r_y = t - t2

    r_out = type(r)((r_x, r_y, r_z), check=False)

    t = p.y + r_z
    t = t ** 2
    t = t - r2
    t = t - r_z ** 2

    t2 = L1 * p.x
    t2 = t2 * 2
    a = t2 - t

    c = r_z * int(q.y * 2)

    b = -L1
    b = b * int(q.x * 2)

    return (a, b, c, r_out)


def line_func_double(r, q):
    """ Requires twisted BN256 point r and BN256 point q. """

    # cache this?
    r_t = r.z ** 2

    A = r.x ** 2
    B = r.y ** 2
    C = B ** 2

    D = r.x + B
    D = D ** 2
    D -= A
    D -= C
    D = D * 2

    E = A * 2 + A
    F = E ** 2

    C8 = C * 8  # C*8

    r_x = F - D * 2
    r_y = E * (D - r_x) - C8

    # (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
    r_z = (r.y + r.z) ** 2 - B - r_t

    # assert r_z == r.y * r.z * 2 # TS/TODO: commented out to test secure computation

    r_out = type(r)((r_x, r_y, r_z), check=False)
    # assert r_out.is_on_curve()

    a = r.x + E
    a = a ** 2
    a -= A + F + B * 4

    t = E * r_t
    t = t * 2
    b = -t
    # print(type(b))
    # print(type(q.x))
    b = b * int(q.x)

    c = r_z * r_t

    c = (c * 2) * int(q.y)

    return (a, b, c, r_out)


def mul_line(r, a, b, c):
    assert type(r) == GFp_12
    assert type(a) == GFp_2
    assert type(b) == GFp_2
    assert type(c) == GFp_2

    # See function fp12e_mul_line in dclxvi

    t1 = GFp_6(gfp_2_zero, a, b)
    t2 = GFp_6(gfp_2_zero, a, b + c)

    t1 = t1 * r.x
    t3 = r.y.mul_scalar(c)
    r.x += r.y
    r.y = t3
    r.x *= t2
    r.x -= t1
    r.x -= r.y
    r.y += t1.mul_tau()


def miller(q, p):
    """ Apply Miller loop on the twisted BN256 point q and BN256 point p.

    Args:
        q, p (CurveElement): twisted and regular BN256 point.

    Returns:
        type(GFp_12): element of tower field GFp_12

    Implementation follows [BDM+10] and [NNS10]:
        [BDM+10]: https://link.springer.com/chapter/10.1007/978-3-642-17455-1_2
        [NNS10]: https://link.springer.com/chapter/10.1007/978-3-642-14712-8_7
    """
    Q = copy.deepcopy(q)
    P = copy.deepcopy(p)
    # TODO: Are Jacobian or affine coordinates required?
    # TODO: Add `assert p.z == 1`? And same for q?

    mQ = ~Q

    f = GFp_12(gfp_6_zero, gfp_6_one)
    T = Q

    Qp = Q.y ** 2

    for naf_i in naf_6up2:
        # Skip on first iteration?
        f = f.square()

        a, b, c, T = line_func_double(T, P)
        mul_line(f, a, b, c)

        if naf_i == 1:
            a, b, c, T = line_func_add(T, Q, P, Qp)
            mul_line(f, a, b, c)
        elif naf_i == -1:
            a, b, c, T = line_func_add(T, mQ, P, Qp)
            mul_line(f, a, b, c)

    Q1 = type(q)((conjugate(Q.x) * xi1[1], conjugate(Q.y) * xi1[2], gfp_2_one))

    Q2 = type(q)((Q.x * (xi2[1].value.value[0]), Q.y, gfp_2_one))

    Qp = Q1.y ** 2
    a, b, c, T = line_func_add(T, Q1, P, Qp)
    mul_line(f, a, b, c)

    Qp = Q2.y ** 2
    a, b, c, T = line_func_add(T, Q2, P, Qp)
    mul_line(f, a, b, c)

    return f


def final_exp(inp):
    assert type(inp) == GFp_12

    # Algorithm 31 from https://eprint.iacr.org/2010/354.pdf

    t1 = inp.conjugate_of()
    inv = inp.inverse()

    t1 = t1.mul(inv)
    # Now t1 = inp^(p**6-1)

    t2 = t1.frobenius_p2()
    t1 = t1.mul(t2)

    fp1 = t1.frobenius()
    fp2 = t1.frobenius_p2()
    fp3 = fp2.frobenius()

    fu1 = t1.exp(u)
    fu2 = fu1.exp(u)
    fu3 = fu2.exp(u)

    y3 = fu1.frobenius()
    fu2p = fu2.frobenius()
    fu3p = fu3.frobenius()
    y2 = fu2.frobenius_p2()

    y0 = fp1.mul(fp2)
    y0 = y0.mul(fp3)

    y1 = t1.conjugate_of()
    y5 = fu2.conjugate_of()
    y3 = y3.conjugate_of()
    y4 = fu1.mul(fu2p)
    y4 = y4.conjugate_of()

    y6 = fu3.mul(fu3p)
    y6 = y6.conjugate_of()

    t0 = y6.square()
    t0 = t0.mul(y4)
    t0 = t0.mul(y5)

    t1 = y3.mul(y5)
    t1 = t1.mul(t0)
    t0 = t0.mul(y2)
    t1 = t1.square()
    t1 = t1.mul(t0)
    t1 = t1.square()
    t0 = t1.mul(y1)
    t1 = t1.mul(y0)
    t0 = t0.square()
    t0 = t0.mul(t1)

    return t0


def optimal_ate(a, b):
    """Optimal ate pairing.

    Original code by Jack Lloyd, see: https://github.com/randombit/pairings.py
    Adapted to integrate with secure_groups and mpyc APIs.

    Invariants:
        a, b are in short weierstrass form.

    Args:
        a (Curve(BN256_TWIST)): BN256_TWIST point
        b (Curve(BN256)): BN256 point

    Returns:
        GFp_12: pairing operation on a, b
    """

    # Convert to Jacobian coordinates with Z==1.
    # NB: a and b can actually be projective. After normalization rest of code will
    # continue as if points are Jacobian.
    a = a.normalize()
    b = b.normalize()

    # if int(a.z) == 0 or int(b.z) == 0:
    #     return GFp_12(gfp_6_zero, gfp_6_one)
    if a is a.identity or b is b.identity:
        return GFp_12(gfp_6_zero, gfp_6_one)

    e = miller(a, b)
    ret = final_exp(e)

    return ret


# TODO:
# 1: Check if conjugate() works if input is 0 or of the form 0i + a; fix
