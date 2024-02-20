""" Implementation of https://eprint.iacr.org/2020/152
``Compressed Σ-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

Protocols:
* Protocol 2, page 13: Pi_s ("pivot")

"""

import os
import sys
import hashlib
from random import SystemRandom
from mpyc.finfields import FiniteFieldElement
import mpyc.mpctools as mpctools
from mpyc.runtime import logging
from mpyc.sectypes import SecureObject
from mpyc.fingroups import EllipticCurvePoint as EllipticCurveElement


prng = SystemRandom()

logger_piv = logging.getLogger("pivot")
logger_piv.setLevel(logging.INFO)

def list_mul(x):
    rettype = type(x[0])
    return mpctools.reduce(rettype.operation, x, initial=rettype.identity)


class AffineForm:
    def __init__(self, coeffs, constant):
        self.coeffs = coeffs
        self.constant = constant

    def __add__(self, other):
        if isinstance(other, AffineForm):
            assert len(self) == len(
                other
            ), "Length of linear forms to add not consistent."
            new_coeffs = [self.coeffs[i] + other.coeffs[i] for i in range(len(self))]
            new_constant = self.constant + other.constant
        elif isinstance(other, (int, FiniteFieldElement, SecureObject)):
            new_coeffs = self.coeffs
            new_constant = self.constant + other
        else:
            raise NotImplementedError(
                f"Addition of form not defined for type: {type(other)}"
            )

        return type(self)(new_coeffs, new_constant)

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __sub__(self, other):
        return self + (-1) * other

    def __mul__(self, other):
        if isinstance(other, (int, FiniteFieldElement)):
            new_coeffs = [coeffs_i * other for coeffs_i in self.coeffs]
            new_constant = self.constant * other
        else:
            raise NotImplementedError(
                f"Multiplication of form not defined for type: {type(other)}"
            )
        return type(self)(new_coeffs, new_constant)

    def __rmul__(self, other):
        return self * other

    def __len__(self):
        return len(self.coeffs)

    def __eq__(self, other):
        return self.coeffs == other.coeffs

    def __repr__(self):
        return f"{str(self.coeffs)}, {str(self.constant)}"

    def eval(self, values):
        assert len(values) == len(
            self.coeffs
        ), "Length of inputs to be equal to coefficients of linear form."
        result = (
            sum([self.coeffs[i] * values_i for i, values_i in enumerate(values)])
            + self.constant
        )
        return result

    def __call__(self, values):
        return self.eval(values)


class LinearForm(AffineForm):
    def __init__(self, coeffs, constant=0):
        self.coeffs = coeffs
        self.constant = 0

    def __add__(self, other):
        if isinstance(other, AffineForm):
            assert len(self) == len(
                other
            ), "Length of linear forms to add not consistent."
            new_coeffs = [self.coeffs[i] + other.coeffs[i] for i in range(len(self))]
            new_constant = self.constant + other.constant
        elif isinstance(other, (int, FiniteFieldElement, SecureObject)):
            new_coeffs = self.coeffs
            new_constant = self.constant + other
        else:
            raise NotImplementedError

        return AffineForm(new_coeffs, new_constant)


def _int(value):
    """ Convert FiniteFieldElements to ints, for sectypes (type: SecureObject) do nothing.

    """
    if isinstance(value, (int, SecureObject)):
        return value
    elif isinstance(value, FiniteFieldElement):
        return int(value)
    else:
        raise NotImplementedError


def fiat_shamir_hash(input_list, order):
    # TODO: Review / Ensure that hash lengths are not ambiguous for group elements.
    # TODO: add to_bytes as a method for Group elements (adds padding to e.g. 200 bytes, not necessarily reversible)
    hash_input = str(input_list).encode("utf-8")
    c = int.from_bytes(hashlib.sha256(hash_input).digest(), "little") % order
    return c


def vector_commitment(x, gamma, g, h):
    """ Pedersen vector commitment. Definition 1 of AC20.
    """
    assert len(g) >= len(x), "Not enough generators."
    prod = list_mul([g[i] ** _int(x_i) for i, x_i in enumerate(x)])
    c = (h ** gamma) * prod
    return c


def affine_to_linear(L, y, n):
    zeros = [0] * n
    constant = L(zeros)
    L_linear = L - constant
    y_linear = y - constant
    return L_linear, y_linear


def prove_linear_form_eval(g, h, P, L, y, x, gamma, gf):
    """ Sigma protocol Pi_s (protocol 2) from AC20.

    Non-interactive version.
    """
    n = len(x)
    L, y = affine_to_linear(L, y, n)
    r = list(gf(prng.randrange(gf.order)) for i in range(n))
    rho = prng.randrange(gf.order)
    t = L(r)
    A = vector_commitment(r, rho, g, h)
    logger_piv.debug(f"Prover computed A={A}.")

    if isinstance(A, EllipticCurveElement):
        input_list = [t, A.normalize(), g, h, P.normalize(), L, y]
    else:
        input_list = [t, A, g, h, P, L, y]

    c = fiat_shamir_hash(input_list, gf.order)
    z = [c * x_i + r[i] for i, x_i in enumerate(x)]
    phi = (c * gamma + rho) % gf.order
    return (
        z,
        phi,
        c,
    )  # TODO: check if it's correct to return c as well (Required to reconstruct A in non-interactive proof)


def verify_linear_form_proof(g, h, P, L, y, z, phi, c):
    n = len(z)
    L, y = affine_to_linear(L, y, n)
    A_check = vector_commitment(z, phi, g, h) * ((P ** c) ** (-1))
    t_check = L(z) - c * y
    logger_piv.debug(f"Verifier computed A_check={A_check}.")
    logger_piv.debug(f"Verifier computed t_check={t_check}.")
    order = type(t_check).order

    if isinstance(A_check, EllipticCurveElement):
        input_list = [t_check, A_check.normalize(), g, h, P.normalize(), L, y]
    else:
        input_list = [t_check, A_check, g, h, P, L, y]

    logger_piv.debug(f"Method verify_linear_form_proof: input_list={input_list}.")
    hash_check = fiat_shamir_hash(input_list, order)
    logger_piv.debug(f"Value of c         ={c}")
    logger_piv.debug(f"Value of hash_check={hash_check}")
    if c == hash_check:
        return True
    else:
        return False
