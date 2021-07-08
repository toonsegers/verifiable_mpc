""" MPC interface for circuit_sat and compressed_pivot modules.

Circuit_sat implementation of https://eprint.iacr.org/2020/152
``Compressed Sigma-Protocol Theory and Practical Application
to Plug & Play Secure Algorithmics''

"""

import mpyc.mpctools as mpctools
import verifiable_mpc.ac20.pivot as pivot
import verifiable_mpc.ac20.circuit_builder as cb
import verifiable_mpc.ac20.circuit_sat_r1cs as cs
from mpyc.runtime import mpc, logging
from verifiable_mpc.ac20.mpc_ac20 import (
    list_mul,
    vector_commitment,
    create_generators,
    koe_trusted_setup,
    koe_restriction_argument_prover,
    koe_opening_linear_form_prover,
    protocol_4_prover,
    protocol_5_prover,
    calculate_fgh_polys,
    _recombination_vectors,
    recombine,
    prove_linear_form_eval,
)

logger_cs_mpc_cb = logging.getLogger("cs_mpc")
logger_cs_mpc_cb.setLevel(logging.INFO)

logger_cs_mpc_cb_hin = logging.getLogger("cs_mpc_hash_inputs")
logger_cs_mpc_cb_hin.setLevel(logging.INFO)

logger_cs_mpc_cb_hout = logging.getLogger("cs_mpc_hash_outputs")
logger_cs_mpc_cb_hout.setLevel(logging.INFO)


async def protocol_8_excl_pivot_prover(generators, circuit, x, gf, use_koe=False):
    """Non-interactive implementation of Protocol 8,
    including nullity protocol, but excluding the call
    to the compressed pivot protocol.

    """
    secfld = type(x[0])
    if "g" in generators:
        g = generators["g"]
        h = generators["h"]
    elif "pp_lhs" in generators:
        use_koe = True
        pp = generators
    else:
        raise NotImplementedError

    logger_cs_mpc_cb.debug(f"Calculate witness.")
    n = len(x)
    assert n == circuit.input_ct
    proof = {}
    m = circuit.mul_ct

    # Calculate a, b, c vectors, only for mul gates
    logger_cs_mpc_cb.debug(f"Calculate a, b, c vectors.")
    a, b, c = circuit.multiplication_triples(x)

    logger_cs_mpc_cb.debug(f"Calculate z.")
    f0 = mpc._random(secfld)
    g0 = mpc._random(secfld)
    a = [f0] + a
    b = [g0] + b
    logger_cs_mpc_cb.debug(f"Calculate z: a, b = await mpc.gather(a,b)")
    a, b = await mpc.gather(a, b)
    logger_cs_mpc_cb.debug(f"Calculate z: fs = recombine(...)")
    fs = recombine(gf, list(zip(range(m + 1), a)), list(range(m + 1, 2 * m + 1)))
    logger_cs_mpc_cb.debug(f"Calculate z: gs = recombine(...)")
    gs = recombine(gf, list(zip(range(m + 1), b)), list(range(m + 1, 2 * m + 1)))
    logger_cs_mpc_cb.debug(f"Calculate z: hs = list(...)")
    hs = list(map(secfld, await mpc.schur_prod(fs, gs)))
    z = x + [f0, g0, f0 * g0] + c + hs

    gamma = mpc._random(secfld)

    if use_koe:
        S = range(len(z))
        z_commitment_P, z_commitment_pi = await koe_restriction_argument_prover(
            S, z, gamma, pp
        )
        z_commitment = {"P": z_commitment_P, "pi": z_commitment_pi}
        proof["z_commitment"] = z_commitment
    else:
        logger_cs_mpc_cb.debug("Calculate commitment for z, denoted by [z].")
        z_commitment = await vector_commitment(z, gamma, g, h)
        # logger_cs_mpc_cb.debug("Calculated [z] in secret shared form.")
        logger_cs_mpc_cb.debug(f"Opened [z] (z_commitment).")
        proof["z_commitment"] = z_commitment

    logger_cs_mpc_cb.debug("Apply first hash of circuit satisfiability protocol.")
    input_list = [z_commitment, str(circuit), "First hash circuit satisfiability protocol"]
    logger_cs_mpc_cb_hin.debug(
        f"Method protocol_8_excl_pivot_prover (1): input_list={input_list}"
    )
    c = pivot.fiat_shamir_hash(
        input_list, gf.order
    )  # TODO: name clash c = a * b from above
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{c}")

    # Prover's response
    linform_f = cb.calculate_fg_form(circuit, wire=0, challenge=c, gf=gf)
    linform_g = cb.calculate_fg_form(circuit, wire=1, challenge=c, gf=gf)
    linform_h = cb.calculate_h_form(circuit, c, gf)

    y1 = linform_f(z)
    y2 = linform_g(z)
    y3 = linform_h(z)
    y1 = await mpc.output(y1, raw = True)  # raw = True ensures that we return field elements
    y2 = await mpc.output(y2, raw = True)
    y3 = await mpc.output(y3, raw = True)

    assert y1 * y2 == y3
    proof["y1"] = y1
    proof["y2"] = y2
    proof["y3"] = y3

    circuit_forms = cb.calculate_circuit_forms(circuit)
    circuit_forms = [cb.convert_to_ac20(f, circuit) for f in circuit_forms]
    outputs = circuit(x)
    outputs = await mpc.output(outputs)
    proof["outputs"] = outputs

    lin_forms = [form - y for form, y in zip(circuit_forms, outputs)] + [
        linform_f - y1,
        linform_g - y2,
        linform_h - y3,
    ]

    # Follow prove_nullity_compressed but add more inputs to Fiat-Shamir hash
    logger_cs_mpc_cb.debug("Apply second hash of circuit satisfiability protocol.")
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
    logger_cs_mpc_cb_hin.debug(
        f"Method protocol_8_excl_pivot_prover (2): input_list={input_list}"
    )
    rho = pivot.fiat_shamir_hash(input_list, gf.order)
    logger_cs_mpc_cb_hout.debug(f"After hash, hash=\n{rho}")
    L = sum((linform_i) * (rho ** i) for i, linform_i in enumerate(lin_forms))
    proof["L"] = L
    return proof, z_commitment, L, z, gamma


async def circuit_sat_prover(
    generators, circuit, x, gf, pivot_choice=cs.PivotChoice.compressed
):
    """Non-interactive implementation of Protocol 8, prover-side,
    including Nullity using compressed pivot (Protocol 5).
    """
    logger_cs_mpc_cb.debug(f"Enter circuit_sat_prover. pivot_choice={pivot_choice}")

    logger_cs_mpc_cb.debug(f"Start protocol 8, excluding pivot proof.")
    proof, z_commitment, L, z, gamma = await protocol_8_excl_pivot_prover(
        generators, circuit, x, gf
    )

    if pivot_choice == cs.PivotChoice.compressed:
        pivot_proof = await protocol_5_prover(
            generators, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == cs.PivotChoice.pivot:
        g = generators["g"]
        h = generators["h"]
        pivot_proof = await prove_linear_form_eval(
            g, h, z_commitment, L, L(z), z, gamma, gf
        )
    elif pivot_choice == cs.PivotChoice.koe:
        L = proof["L"]
        P = proof["z_commitment"]["P"]
        pi = proof["z_commitment"]["pi"]
        pivot_proof, u = await koe_opening_linear_form_prover(L, z, gamma, generators, P, pi)
    else:
        raise NotImplementedError
    proof["pivot_proof"] = pivot_proof

    return proof


