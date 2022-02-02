import verifiable_mpc.tools.code_to_r1cs as c2r
import verifiable_mpc.tools.qap_creator as qc


class QAP:

    def __init__(self, code, field):
        inputs, body = c2r.extract_inputs_and_body(c2r.parse(code))
        flatcode = c2r.flatten_body(body)
        varnames = c2r.get_var_placement(inputs, flatcode)
        V, W, Y = c2r.flatcode_to_r1cs(inputs, flatcode)
        V = for_each_in(int, lambda x: field(x), V)
        W = for_each_in(int, lambda x: field(x), W)
        Y = for_each_in(int, lambda x: field(x), Y)
        v, w, y, t = qc.r1cs_to_qap_ff(V, W, Y, field)
        v = [qc.Poly(coeffs) for coeffs in v]
        w = [qc.Poly(coeffs) for coeffs in w]
        y = [qc.Poly(coeffs) for coeffs in y]
        t = qc.Poly(t)
        self.v = v
        self.w = w
        self.y = y
        self.t = t
        self.inputs = inputs
        self.flatcode = flatcode
        self.varnames = varnames
        self.d = len(flatcode)
        self.m = len(varnames) - 1 # Note: `~one` is not included in count
        self.out_ix = varnames.index("~out")
        # TODO: allow for multiple output indices (previous line, also below)
        self.indices = range(self.m+1) 
        self.indices_io_and_0 = range(0, self.out_ix + 1) # includes "one" 
        self.indices_io = range(1, self.out_ix + 1) # io indices exclude "one" 
        self.indices_mid = range(self.out_ix+1, self.m+1)

    def calculate_witness(self, input_vars):
        witness = c2r.assign_variables(self.inputs, input_vars, self.flatcode)
        assert int(witness[0]) == 1, "First coordinate of witness != 1"
        return witness


# Helper method from pysnark.runtime (credits to Meilof Veeningen)
def for_each_in(cls, f, struct):
    """ Recursively traversing all lists and tuples in struct, apply f to each
        element that is an instance of cls. Returns structure with f applied. """
    if isinstance(struct, list):
        return list(map(lambda x: for_each_in(cls, f, x), struct))
    elif isinstance(struct, tuple):
        return tuple(map(lambda x: for_each_in(cls, f, x), struct))
    else:
        if isinstance(struct, cls):
            return f(struct)
        else:
            return struct


def calculate_witness(code, input_vars):
    inputs, body = c2r.extract_inputs_and_body(c2r.parse(code))
    flatcode = c2r.flatten_body(body)
    witness = c2r.assign_variables(inputs, input_vars, flatcode)
    return witness
