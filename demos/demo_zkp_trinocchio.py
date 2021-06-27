"""Demo Trinocchio protocol based on pynocchio module.

Credits:
* Trinocchio by Berry Schoenmakers, Meilof Veeningen and 
  Niels de Vreede: https://eprint.iacr.org/2015/480
* MPyC by Berry Schoenmakers: https://github.com/lschoe/mpyc
* bn256 by Jack Lloyd: https://github.com/randombit/pairings.py/ 
  (BSD-2-Clause license)
* r1cs and qap tools by Vitalik Buterin: 
  https://github.com/ethereum/research/tree/master/zksnark (MIT license)
"""

import os, sys
import pprint as pp

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import verifiable_mpc.trinocchio.pynocchio as pynocchio
import verifiable_mpc.trinocchio.trinocchio as trinocchio
from mpyc.runtime import mpc
from mpyc.thresha import _recombination_vector
import verifiable_mpc.tools.code_to_r1cs_and_qap.code_to_qap as c2q
import verifiable_mpc.tools.code_to_r1cs_and_qap.qap_creator as qc


async def main():
    await mpc.start()
    print(f"MPC parties {mpc.parties}")
    m = len(mpc.parties)
    trusted_party_id = 0
    print(f"Trusted party PID: {trusted_party_id}")

    secfld = mpc.SecFld(modulus=trinocchio.modulus)
    gf = secfld.field
    gf.is_signed = False

    # Client provides function code and inputs

    #     inputs = [secfld(3)]
    #     code = """
    # def qeval(x):
    #     y = x**3
    #     return y + x + 5
    # """
    inputs = [secfld(3), secfld(2)]
    code = """
def qeval(x, y):
    z = x**3 + 2*y**2
    return z + x + 5
"""

    # QAP creation step
    qap = c2q.QAP(code, gf)
    print(f"QAP created. Size: {qap.m}, degree {qap.d}.")

    # Trusted party's KeyGen step
    if mpc.pid == trusted_party_id:
        td = pynocchio.Trapdoor(trinocchio.modulus)
        gen = pynocchio.Generators(td, trinocchio.g1, trinocchio.g2)
        evalkey = pynocchio.generate_evalkey(td, qap, gen)
        verikey = pynocchio.generate_verikey(td, qap, gen)
    else:
        evalkey = None
        verikey = None
    print("Keysets generated.")

    # Trusted party communicates keysets
    evalkey = await mpc.transfer(evalkey, trusted_party_id)
    verikey = await mpc.transfer(verikey, trusted_party_id)
    print("Trusted setup completed. Keysets received by parties.")

    # Prover's steps
    c = qap.calculate_witness(inputs)
    p = pynocchio.compute_p_poly(qap, c)
    h, r = p / qap.t

    # Compute proof_input for this MPC party
    # TODO: Proof not zero knowledge (not required for verification to work); make it ZK
    c_shares = await mpc.gather(c)
    # Compute shares of coefficients of h polynomial
    h_coeffs_shares = await mpc.gather(h.coeffs)
    h_shares = qc.Poly(h_coeffs_shares)
    proof_input = pynocchio.compute_proof(qap, c_shares, h_shares, evalkey)
    print("Proof computed.")

    # Communicate proof elements to all parties
    proof_inputs = await mpc.transfer(proof_input)

    # Recombine proof elements
    xcoords = tuple(i + 1 for i in range(m))
    lagrange_vect = _recombination_vector(gf, xcoords, 0)
    proof = {}
    for key in proof_input.keys():
        points_lambda = [proof_inputs[i][key] * lagrange_vect[i] for i in range(m)]
        proof_element = pynocchio.apply_to_list(trinocchio.point_add, points_lambda)
        proof[key] = proof_element
    print("Proof recombined.")

    c_out = await mpc.output(c[1:])
    c_out = [1] + c_out
    c_client = c_out[: qap.out_ix + 1]

    verifications = pynocchio.verify(qap, verikey, proof, c_client)
    if all(verifications.values()):
        print("All checks passed.")
    else:
        print("Not all checks passed.")
    pp.pprint(verifications)

    await mpc.shutdown()


if __name__ == "__main__":
    mpc.run(main())


# TODO list
# TODO-1: Implement ZK variant
# TODO-1: Implement multi-client with basic Trinocchio
# TODO-1: Implement full Trinocchio scheme (Alg. 4, p. 25)
