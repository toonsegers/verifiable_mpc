"""Implements https://eprint.iacr.org/2020/152.

This module implements the circuit satisfiability protocol
(Protocol 8) from ``Compressed Sigma-Protocol Theory and
Practical Application to Plug & Play Secure Algorithmics''.

Example:
    >>> proof = circuit_sat_prover(generators, circuit, x, ...)
    >>> tests = circuit_sat_verifier(proof, ...)

"""

import os
import sys
# from enum import Enum
import logging
from random import SystemRandom
# import re

import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.compressed_pivot as compressed_pivot
import verifiable_mpc.ac20.knowledge_of_exponent as koe
from verifiable_mpc.ac20.circuit_sat_r1cs import (
    PivotChoice,
    create_generators,
    next_power_of_2,
    lagrange,
    calculate_fgh_polys,
)

# from mpyc.thresha import _recombination_vector as lagrange   ## see below
import verifiable_mpc.ac20.circuit_builder as cb

prng = SystemRandom()

logger_cs2 = logging.getLogger("circuit_sat")
logger_cs2.setLevel(logging.INFO)

logger_cs2_hin = logging.getLogger("circuit_sat_hash_inputs")
logger_cs2_hin.setLevel(logging.INFO)

logger_cs2_hout = logging.getLogger("circuit_sat_hash_outputs")
logger_cs2_hout.setLevel(logging.INFO)


def check_input_length_power_of_2(x, circuit, padding_value=0):
    """Pad code and x to meet requirements of compressed pivot protocol."""
    assert circuit.input_ct == len(x)
    z_len = circuit.input_ct + 3 + 2 * circuit.mul_ct
    # Calculate required padding of x such that len(z) + 1 is a power of 2
    if not bin(z_len + 1).count("1") == 1:
        padding = next_power_of_2(z_len) - z_len - 1
    else:
        padding = 0
    check = padding == 0
    return check, padding, z_len + padding


def protocol_8_excl_pivot_prover(generators, circuit, x, gf, use_koe=False):
    """Non-interactive implementation of Protocol 8,
    including nullity protocol, but excluding the call
    to the compressed pivot protocol.

    """
    if "g" in generators:
        g = generators["g"]
        h = generators["h"]
    elif "pp_lhs" in generators:
        use_koe = True
        pp = generators
    else:
        raise NotImplementedError

    n = len(x)
    assert n == circuit.input_ct
    proof = {}

    m = circuit.mul_ct

    # Calculate a, b, c vectors, only for mul gates
    a, b, c = circuit.multiplication_triples(x)

    # Calculate random polynomials f, g, h
    f_poly, g_poly, h_poly = calculate_fgh_polys(
        a, b, None, gf
    )  # TODO: remove third param from function

    # Prover commits to (x, f(0), g(0), h(0), h(1), .., h(2m))
    h_evaluations = [h_poly.eval(i + 1) for i in range(2 * m)]
    z = x + [f_poly.eval(0), g_poly.eval(0), h_poly.eval(0)] + h_evaluations

    gamma = prng.randrange(1, gf.order)

    if use_koe:
        S = range(len(z))
        z_commitment_P, z_commitment_pi = koe.restriction_argument_prover(
            S, z, gamma, pp
        )
        z_commitment = {"P": z_commitment_P, "pi": z_commitment_pi}
        proof["z_commitment"] = z_commitment
    else:
        logger_cs2.debug("Calculate [Z].")
        z_commitment = pivot.vector_commitment(z, gamma, g, h)
        proof["z_commitment"] = z_commitment

    logger_cs2.debug("Apply first hash of circuit satisfiability protocol.")
    input_list = [z_commitment, str(circuit), "First hash circuit satisfiability protocol"]
    logger_cs2_hin.debug(
        f"Method protocol_8_excl_pivot_prover: input_list=\n{input_list}"
    )
    c = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs2_hout.debug(f"After hash, hash=\n{c}")

    # Prover's response
    y1 = f_poly.eval(c)
    y2 = g_poly.eval(c)
    y3 = h_poly.eval(c)
    assert y3 == y1 * y2

    """ Create linear forms corresponding to f(c), g(c), h(c) for Pi_Nullity step
    in Protocol 6 (Protocol 8 in updated version)
    """
    linform_f = cb.calculate_fg_form(circuit, wire=0, challenge=c, gf=gf)
    linform_g = cb.calculate_fg_form(circuit, wire=1, challenge=c, gf=gf)
    linform_h = cb.calculate_h_form(circuit, c, gf)

    y1 = linform_f(z)
    y2 = linform_g(z)
    y3 = linform_h(z)
    assert y1 * y2 == y3
    proof["y1"] = y1
    proof["y2"] = y2
    proof["y3"] = y3

    # Calculate linear form(s) for circuit, and output value(s)
    circuit_forms = cb.calculate_circuit_forms(circuit)
    circuit_forms = [cb.convert_to_ac20(f, circuit) for f in circuit_forms]
    outputs = circuit(x)
    proof["outputs"] = outputs

    lin_forms = [form - y for form, y in zip(circuit_forms, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    # Follow prove_nullity_compressed but add more inputs to Fiat-Shamir hash
    logger_cs2.debug("Apply second hash of circuit satisfiability protocol.")
    input_list = [
        y1,
        y2,
        y3,
        z_commitment,
        outputs,
        circuit_forms,
        lin_forms,
        "Second hash circuit satisfiability protocol",
    ]
    logger_cs2_hin.debug(
        f"Method: protocol_8_excl_pivot_prover, input_list=\n{input_list}"
    )
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs2_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    proof["L"] = L
    return proof, z_commitment, L, z, gamma


def protocol_8_excl_pivot_verifier(proof, circuit, gf, use_koe=False):
    verification = {}
    y1 = proof["y1"]
    y2 = proof["y2"]
    y3 = proof["y3"]
    if not y1 * y2 == y3:
        verification["y1*y2=y3"] = False
        return verification
    else:
        verification["y1*y2=y3"] = True

    n = circuit.input_ct
    m = circuit.output_ct

    if "P" in proof:
        use_koe = True
    if use_koe:
        z_commitment_P = proof["z_commitment"]["P"]
        z_commitment_pi = proof["z_commitment"]["pi"]
        logger_cs2.debug("Apply first hash of circuit satisfiability protocol.")
        input_list = [
            z_commitment_P,
            z_commitment_pi,
            str(circuit),
            "First hash circuit satisfiability protocol",
        ]
        logger_cs2_hin.debug(
            f"Method: protocol_8_excl_pivot_verifier (1), input_list=\n{input_list}"
        )
        c = pivot.fiat_shamir_hash(input_list, gf.order)
        logger_cs2_hout.debug(f"After hash, hash=\n{c}")
    else:
        z_commitment = proof["z_commitment"]
        logger_cs2.debug("Apply first hash of circuit satisfiability protocol.")
        input_list = [
            z_commitment,
            str(circuit),
            "First hash circuit satisfiability protocol",
        ]
        logger_cs2_hin.debug(
            f"Method: protocol_8_excl_pivot_verifier (1), input_list=\n{input_list}"
        )
        c = pivot.fiat_shamir_hash(input_list, gf.order)
        logger_cs2_hout.debug(f"After hash, hash=\n{c}")

    linform_f = cb.calculate_fg_form(circuit, wire=0, challenge=c, gf=gf)
    linform_g = cb.calculate_fg_form(circuit, wire=1, challenge=c, gf=gf)
    linform_h = cb.calculate_h_form(circuit, c, gf)

    outputs = proof["outputs"]
    circuit_forms = cb.calculate_circuit_forms(circuit)
    circuit_forms = [cb.convert_to_ac20(f, circuit) for f in circuit_forms]

    lin_forms = [form - output for form, output in zip(circuit_forms, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    logger_cs2.debug("Apply second hash of circuit satisfiability protocol.")
    input_list = [
        y1,
        y2,
        y3,
        z_commitment,
        outputs,
        circuit_forms,
        lin_forms,
        "Second hash circuit satisfiability protocol",
    ]
    logger_cs2_hin.debug(
        f"Method: protocol_8_excl_pivot_verifier (2), input_list=\n{input_list}"
    )
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs2_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))

    if not L == proof["L"]:
        verification["L_wellformed_from_Cfgh_forms"] = False
        return verification, L
    else:
        verification["L_wellformed_from_Cfgh_forms"] = True
    
    return verification, L


def circuit_sat_prover(generators, circuit, x, gf, pivot_choice=PivotChoice.compressed):
    """Non-interactive implementation of Protocol 8, prover-side,
    including Nullity using arbitrary pivot.
    """
    proof, z_commitment, L, z, gamma = protocol_8_excl_pivot_prover(
        generators, circuit, x, gf
    )

    if pivot_choice == PivotChoice.compressed:
        pivot_proof = compressed_pivot.protocol_5_prover(
            generators, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == PivotChoice.pivot:
        g = generators["g"]
        h = generators["h"]
        pivot_proof = pivot.prove_linear_form_eval(
            g, h, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == PivotChoice.koe:
        L = proof["L"]
        P = proof["z_commitment"]["P"]
        pi = proof["z_commitment"]["pi"]
        pivot_proof, u = koe.opening_linear_form_prover(L, z, gamma, generators, P, pi)
    else:
        raise NotImplementedError
    proof["pivot_proof"] = pivot_proof

    return proof


def circuit_sat_verifier(
    proof, generators, circuit, gf, pivot_choice=PivotChoice.compressed
):
    """Non-interactive implementation of Protocol 8, verifier-side,
    including Nullity using arbitrary pivot.
    """
    verification, L = protocol_8_excl_pivot_verifier(proof, circuit, gf)

    # Run pivot protocol based on enum `pivot_choice`
    if pivot_choice == PivotChoice.compressed:
        z_commitment = proof["z_commitment"]
        pivot_proof = proof["pivot_proof"]
        pivot_verification = compressed_pivot.protocol_5_verifier(
            generators, z_commitment, L, 0, pivot_proof, gf
        )
    elif pivot_choice == PivotChoice.pivot:
        z_commitment = proof["z_commitment"]
        g = generators["g"]
        h = generators["h"]
        z, phi, c = proof["pivot_proof"]
        pivot_verification = pivot.verify_linear_form_proof(
            g, h, z_commitment, L, 0, z, phi, c
        )
    elif pivot_choice == PivotChoice.koe:
        koe_pivot_proof = proof["pivot_proof"]
        pivot_verification = koe.opening_linear_form_verifier(
            L, generators, koe_pivot_proof, 0
        )
        verification["pivot_verification"] = pivot_verification
    else:
        raise NotImplementedError
    verification["pivot_verification"] = pivot_verification

    return verification


# TODO-list
# TODO-1: Avoid double work: compute z_commitment once, don't compute P (equal to z_commitment) again
# TODO-1: Apply DRY to circuit_sat sub-methods; remove redundant code
# TODO-1: Rename `symbols` to `vars`
# TODO-1: Clean up create_fgh_linear_forms
# TODO-1: Add "tag" to hashes
# TODO-3: Reduce mul-gates: represent multiply-by-constant lines in flatcode by constant in R1CS matrices
