# WORK-IN-PROGRESS

""" Implementation of Geppetri key gen steps with BN256 curves 

Using pairings.py from https://github.com/randombit/pairings.py
Alternative for bn128 by Ethereum: https://github.com/ethereum/py_ecc/blob/master/py_ecc/bn128/bn128_curve.py

"""

import os, sys
from random import SystemRandom

import zk_helpers.pairings.bn256 as bn256


# k,g = bn256.g1_random()
# k,g = bn256.g2_random()

# in qap2key:
# void generate_master_skey(mastersk& sk) {
#     sk.s = modp::rand();
#     sk.rc = modp::rand();
#     sk.al = modp::rand();
# }

DEFAULT_QAP_DEGREE = 20

def list_add(grp_elements):
    n = len(grp_elements)
    if n == 1:
        return grp_elements[0]

    m0 = list_add(grp_elements[:n//2])  
    m1 = list_add(grp_elements[n//2:])
    return m0.add(m1)


def generate_s():
    """ Generate secret s """
    s = prng.randrange(bn256.order)
    return s


def generate_crs(s, qap_degree = DEFAULT_QAP_DEGREE):
    """ Generate common reference string per function G01, p. 7 """
    crs_g1 = {"x^"+str(i)+"_g1"  : g1.scalar_mul(s**i) for i in range(qap_degree + 1)}
    crs_g2 = {"x^"+str(i)+"_g2"  : g2.scalar_mul(s**i) for i in range(qap_degree + 1)}
    return {**crs_g1, **crs_g2}


def generate_commitment_key(qap_degree = DEFAULT_QAP_DEGREE):
    """ Generate commitment key per function Gc1, p. 7 """
    alpha = prng.randrange(bn256.order)
    ck_g1 = {"x^"+str(i)+"_g1": g1.scalar_mul(s**i) for i in range(qap_degree + 1)}
    ck_g2 = {"ax^"+str(i)+"_g2": g2.scalar_mul(alpha * (s**i)) for i in range(qap_degree + 1)}
    return {**ck_g1, **ck_g2}


def commit(v, r, ck):
    """ Commit to input v per function C1, p. 7 """
    c_g1 = ck["x^"+str(0)+"_g1"].scalar_mul(r) 
    x_terms = list_add([ck["x^"+str(i+1)+"_g1"].scalar_mul(v[i]) for i in range(len(v))])
    c_g1 = c_g1.add(x_terms)

    c_g2 = ck["ax^"+str(0)+"_g2"].scalar_mul(r) 
    x_terms = list_add([ck["ax^"+str(i+1)+"_g2"].scalar_mul(v[i]) for i in range(len(v))])
    c_g2 = c_g2.add(x_terms)

    return (c_g1, c_g2)


def code_to_qap(code, ffield):
    r1cs = code_to_r1cs(code, ffield)
    qap = r1cs_to_qap(r1cs, ffield)
    return qap # v, w, y, t



if __name__ == '__main__':
    prng = SystemRandom()
    g1 = bn256.curve_G
    g2 = bn256.twist_G
    d = 5 # degree of the QAP
    n = 5 # input length

    s = generate_s()
    crs = generate_crs(s, d)
    ck = generate_commitment_key(d)
    v = [1, 2, 3, 4, 5]
    r = prng.randrange(bn256.order)
    c = commit(v, r, ck)
    print(c)
    

