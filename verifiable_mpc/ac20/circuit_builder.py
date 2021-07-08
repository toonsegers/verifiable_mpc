from enum import Enum

from mpyc.finfields import GF, FiniteFieldElement, PrimeFieldElement
from mpyc.runtime import mpc
from mpyc.sectypes import SecureFiniteField, SecureInteger
import verifiable_mpc.ac20.mpc_ac20
from verifiable_mpc.ac20.pivot import AffineForm, LinearForm
from verifiable_mpc.ac20.circuit_sat_r1cs import calculate_fgh_polys
import verifiable_mpc.tools.qap_creator as qc
from random import SystemRandom


class op(Enum):
    add = "add"
    mul = "mul"
    scalar_mul = "scalar mul"


class Gate:
    """Note: AC20 requires fan-in 2 (unbounded fan-out)"""

    def __init__(self, op, output, inputs):
        self.op = op
        # TODO: allow for multiple outputs
        self.output = output  # Includes output variables
        self.inputs = inputs  # Includes input variables
        self.mul_index = None
        self.index = None

    def __str__(self):
        inputs = str([i.name if isinstance(i, CircuitVar) else i for i in self.inputs])
        output = str(
            self.output.name if isinstance(self.output, CircuitVar) else self.output
        )
        return output + " <- " + str(self.op) + "(" + inputs + ")"


class Circuit:
    def __init__(self):
        self.gates = []
        self.gate_ct = 0
        self.input_ct = 0  # Count input variables
        self.output_ct = 0  # Count output variables
        self.add_ct = 0  # Count add gates
        self.mul_ct = 0  # Count mul gates
        self.scalar_mul_ct = 0  # Count scalar_mul gates
        self._dummy_ct = 0  # Count dummy gates
        self.input_gates = []  # Track indices
        self.output_gates = []  # Track indices
        self.circuitvars = []

    def add_gate(self, gate):
        """Add gate to circuit."""
        self.gates.append(gate)
        self.gate_ct += 1
        gate.index = self.gate_ct - 1

        # Add gate.index to input CircuitVars. Add CircuitVar to self.circuitvars. (If not present in list.)
        for i in [0, 1]:
            if isinstance(gate.inputs[i], CircuitVar):
                if gate.index not in gate.inputs[i].gates:
                    gate.inputs[i].gates.append(gate.index)

        # If input-wire CircuitVars are circuit inputs, then append to circuit.input_gates list
        if (
            isinstance(gate.inputs[0], CircuitVar)
            and gate.inputs[0].input_index is not None
        ):
            self.input_gates.append(gate.index)
        elif (
            isinstance(gate.inputs[1], CircuitVar)
            and gate.inputs[1].input_index is not None
        ):
            self.input_gates.append(gate.index)

        # Increment add, mul, scalar_mul counters and give gate a mul_index if needed.
        if gate.op == op.add:
            self.add_ct += 1
        elif gate.op == op.mul:
            assert isinstance(gate.inputs[0], CircuitVar) and isinstance(
                gate.inputs[1], CircuitVar
            )
            self.mul_ct += 1
            # Give each mul-gate an index
            gate.mul_index = self.mul_ct - 1  # index := count-1
        elif gate.op == op.scalar_mul:
            self.scalar_mul_ct += 1
        else:
            raise NotImplementedError

    def name_dummy(self):
        """Give new circuit variable a dummy name."""
        name = "dummy_" + str(self._dummy_ct)
        self._dummy_ct += 1
        return name

    # TODO: change to property
    def parents(self, gate):
        """Returns parent gates of gate."""
        parents = [
            g
            for g in self.gates
            if gate.output.name
            in [v.name for v in g.inputs if isinstance(v, CircuitVar)]
        ]
        return parents

    # TODO: change to property
    def children(self, gate):
        """Returns children gates of gate."""
        gate_inputs = [v.name for v in gate.inputs if isinstance(v, CircuitVar)]
        children = [g for g in self.gates if g.output.name in gate_inputs]
        return children

    # TODO: change to property
    def mul_gates(self):
        return [g for g in self.gates if g.op == op.mul]

    def out_gates(self):
        """Returns output gates (not indices of output gates)."""
        return [self.gates[ix] for ix in self.output_gates]

    def in_gates(self):
        """Returns input gates (not indices of input gates)."""
        return [self.gates[ix] for ix in self.input_gates]

    def initial_inputs(self):
        # TODO: Ensure that list is always ordered by .input_index
        return [v.value for v in self.circuitvars if v.input_index != None]

    def multiplication_triples(self, inputs):
        """Return left, right and output wire values for all mul-gates."""
        left_wire_forms = [
            construct_affine_form(g, self, wire=0)
            for i, g in enumerate(self.mul_gates())
        ]
        right_wire_forms = [
            construct_affine_form(g, self, wire=1)
            for i, g in enumerate(self.mul_gates())
        ]
        wire_forms = zip(left_wire_forms, right_wire_forms)
        alpha = [0] * self.mul_ct
        beta = [0] * self.mul_ct
        gamma = [0] * self.mul_ct
        for i, (left, right) in enumerate(wire_forms):
            alpha[i] = left(inputs + gamma)
            beta[i] = right(inputs + gamma)
            gamma[i] = alpha[i] * beta[i]
        return alpha, beta, gamma

    def eval(self, inputs, gate):
        """Evaluate circuit up to given gate."""
        # inputs = inputs[:-1*self.input_padding]  # If we applied padding to the input vector, remove it.
        _, _, gamma = self.multiplication_triples(inputs)
        form_gate_l = construct_affine_form(gate, self, wire=0)
        form_gate_r = construct_affine_form(gate, self, wire=1)
        left = form_gate_l(inputs + gamma)
        right = form_gate_r(inputs + gamma)
        if gate.op == op.add:
            y = left + right
        elif gate.op in [op.mul, op.scalar_mul]:
            y = left * right
        else:
            raise ValueError
        return y

    def __call__(self, inputs):
        """Evaluate circuit for all output gates."""
        y = [self.eval(inputs, self.gates[gate_ix]) for gate_ix in self.output_gates]
        return y

    def __str__(self):
        return print_circuit(self)


class CircuitVar:
    """Wrap type with methods, attributes that record circuit attributes (variables, gates).

    Inputs:
        value       input instance
        circuit     Circuit instance
        name        name of variable
        input_var   flag to indicate input variable (for circuit.input_ct counter)
        output_var  flag to indicate output variable (for circuit.output_ct counter)

    """

    def __init__(self, value, circuit, name=None, input_var=True):
        self.value = value
        self.circuit = circuit
        self.name = name
        self.input_index = None  # if index != None, then self is an input
        self.output_index = None  # if index != None, then self is an output
        self.gates = []  # Track indices of gates that use given CircuitVar

        # Track input variables, index them
        if input_var:
            circuit.input_ct += 1
            self.input_index = circuit.input_ct - 1  # index := count-1
            self.name += "_input_" + str(self.input_index)

        circuit.circuitvars.append(self)

    def label_output(a, name):
        # Allow user to index variable as output after instantiation.
        if a.output_index == None:
            a.circuit.output_ct += 1
            a.output_index = a.circuit.output_ct - 1  # index := count-1
            if name:
                a.name = name + "_output_" + str(a.output_index)
            else:
                a.name = a.name + "_output_" + str(a.output_index)

        # Add gates that include a to circuit.output_gates
        output_gates = [g.index for g in a.circuit.gates if g.output is a]
        a.circuit.output_gates.extend(output_gates)

    def __add__(self, right):
        if isinstance(right, CircuitVar):
            value = self.value + right.value
        elif isinstance(right, (int, FiniteFieldElement)):
            value = self.value + right
        else:
            raise NotImplementedError
        out = type(self)(
            value, self.circuit, name=self.circuit.name_dummy(), input_var=False
        )
        g = Gate(op.add, out, [self, right])
        self.circuit.add_gate(g)
        return out

    def __radd__(self, right):
        return self + right

    def __sub__(self, right):
        return self + (-1 * right)

    def __rsub__(self, right):
        return (-1 * self) + right

    def __mul__(self, right):
        if isinstance(right, CircuitVar):
            value = self.value * right.value
        elif isinstance(right, (int, FiniteFieldElement)):
            value = self.value * right
        else:
            raise NotImplementedError

        out = type(self)(
            value, self.circuit, name=self.circuit.name_dummy(), input_var=False
        )
        if isinstance(right, CircuitVar):
            g = Gate(op.mul, out, [self, right])
        elif isinstance(right, (int, FiniteFieldElement)):
            g = Gate(op.scalar_mul, out, [self, right])
        else:
            raise NotImplementedError
        self.circuit.add_gate(g)
        return out

    def __rmul__(self, right):
        return self * right

    def check_not_zero(self):
        """Gadget that implements b = (a != 0) ? 1 : 0.

        Computes output, witnesses and required gates.
        """

        # Compute output CircuitVar
        a = self.value
        b = mpc.if_else(a == 0, 0, 1)

        # Calculate witnesses
        assert isinstance(a, (FiniteFieldElement, SecureFiniteField, int, SecureInteger))
        c = (a + (1 - b)) ** (-1)
        cv_c = type(self)(c, self.circuit, name="witness_{" + self.name + "!=0}", input_var=True)

        # Expand circuit with gates for equations a · c = b, a · (1 − b) = 0
        cv_b = self * cv_c
        cv_d = self * (1 - cv_b)
        cv_d.label_output("witness_{" + self.name + "!=0}") 

        return cv_b

    def __ne__(self, other):
        return (self - other).check_not_zero()

    def check_ge_zero(self):
        """Gadget that implements b = (a >= 0) ? 1 : 0.

        Computes output, witnesses and required gates.
        """
        # TODO: test new gate (for secint and int)
        # Compute output CircuitVar
        a = self.value
        b = a >= 0
        cv_b = type(self)(b, self.circuit, name="witness_{" + self.name + ">=0}", input_var=True)

        # Calculate witnesses
        assert isinstance(a, (int, SecureInteger))
        if isinstance(a, SecureInteger):
            c = mpc.to_bits(a)
        elif isinstance(a, int):
            # TODO: review
            c = twos_complement(a, a.bit_length() + 1)
        else:
            raise NotImplementedError

        cv_c = [
            type(self)(
                c_i, self.circuit, name="witness_{" + self.name + ">=0}", input_var=True
            )
            for c_i in c
        ]

        # Expand circuit with gates for equation d = sum c_i * 2^i (in twos complement) and c_i * c_i = c_i
        cv_d = -1 * cv_c[-1] * 2 ** (len(cv_c) - 1) + sum(
            cv_c_i * 2 ** i for i, cv_c_i in enumerate(cv_c[:-1])
        )
        # Prover uses circuit to show that d = sum c_i * 2^i (in twos complement)
        cv_d.label_output("witness_{" + self.name + ">=0}")

        e = [cv_c_i * cv_c_i - cv_c_i for cv_c_i in cv_c]
        [e_i.label_output("witness_{" + self.name + ">=0}") for e_i in e]

        return cv_b

    def __le__(self, other):
        return (other - self).check_ge_zero()

    def __lt__(self, other):
        return (other - self - 1).check_ge_zero()

    def __gt__(self, other):
        return (self - other - 1).check_ge_zero()

    def __ge__(self, other):
        return (self - other).check_ge_zero()

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return self.name + "{" + str(self.value) + "}"


def twos_complement(value, bit_length):
    return bin(value & (2 ** bit_length - 1))


def print_circuit(circuit):
    ret = ""
    for gate in circuit.out_gates():
        ret += print_out_gate(circuit, gate)
    return ret


def print_out_gate(circuit, gate, level=0):
    ret = "\t" * level + str(gate) + "\n"
    for child in circuit.children(gate):
        ret += print_out_gate(circuit, child, level + 1)
    return ret


def construct_affine_form(gate, circuit, wire=None):
    """Construct affine form for left (resp. right) wire of gate.

    Inputs:
        gate        Gate instance indicating which mul-gate
        circuit     Circuit instance indicating circuit to traverse
        wire        0 indicates left, 1 indicates right wire, None indicates all

    Output:
        AffireForm of length circuit.input_ct + circuit.mul_ct
        Note: AC20 requires AffineForm of length input_ct + 3 + 2*mul_ct. Convert later.
    """

    def construct_for_wire(gate, circuit, wire):
        """Update given affine form for left or right input wire."""
        assert wire in [0, 1]  # wire is either left or right

        # Start with zero
        ret = AffineForm([0] * circuit.input_ct + [0] * circuit.mul_ct, 0)

        # If input is a constant (not a Circuit_Var)
        if not isinstance(gate.inputs[wire], CircuitVar):
            ret.constant += gate.inputs[wire]
        # If input is Circuit_Var
        else:
            # If input is circuit input
            if gate.inputs[wire].input_index is not None:
                ret.coeffs[gate.inputs[wire].input_index] += 1
            # Input is mul-gate or add-gate
            else:
                gate_input = gate.inputs[wire].name
                child_gate = [g for g in circuit.gates if g.output.name == gate_input][0]
                # If input is mul-gate
                if child_gate.op == op.mul:
                    ret.coeffs[circuit.input_ct + child_gate.mul_index] += 1
                # If input is add-gate
                elif child_gate.op == op.add:
                    ret = construct_affine_form(child_gate, circuit, wire=None)
                # If input is scalar_mul-gate
                elif child_gate.op == op.scalar_mul:
                    ret = construct_affine_form(child_gate, circuit, wire=None)
                else:
                    # All valid options are tested at this point.
                    raise ValueError
        return ret

    # Start with all-zero affine form, also when called from construct_for_wire
    ret = AffineForm([0] * circuit.input_ct + [0] * circuit.mul_ct, 0)

    # Input is mul-gate and designation to traverse left or right wire
    if wire is not None:
        # Construct affine form for specified input wire (pass all-zero affine form)
        ret = construct_for_wire(gate, circuit, wire)
    # wire == None; we are in the recursion and need to traverse left and right wires, or we have passed an output wire
    elif wire == None:
        if gate.op == op.add:
            for _wire in [0, 1]:
                ret += construct_for_wire(gate, circuit, _wire)
        elif gate.op == op.scalar_mul:
            # left wire is CircuitVar, then right wire should be scalar
            if isinstance(gate.inputs[0], CircuitVar):
                ret = construct_for_wire(gate, circuit, 0)
                ret *= gate.inputs[1]
            # right wire is CircuitVar, then left wire should be scalar
            elif isinstance(gate.inputs[1], CircuitVar):
                ret = construct_for_wire(gate, circuit, 1)
                ret *= gate.inputs[0]
            # no input wires connect to a CircuitVar. Hence, both are scalars
            else:
                s = gate.inputs[0] * gate.inputs[1]
                ret.constant = s
        elif gate.op == op.mul:
            # If we passed an output gate that is a mul-gate, then set coefficient corresponding to mul-gate to 1.
            assert gate.output.output_index != None
            ret.coeffs[circuit.input_ct + gate.mul_index] = 1
            # Mul-gate is not possible in the recursion: construct_for_wire cannot call construct_affine_form
            # for mul-gate, it should return an affine form with increased mul-index
        else:
            raise ValueError
    else:
        raise ValueError
    return ret


def convert_to_ac20(form, circuit):
    """Include 0 coefficients for f(0), g(0), h(0) and h(i), i = m+1, ..., 2m.

    Ensures that input vector to form() corresponds to z-vector in AC20 of length n + 3 + 2m.
    """
    newform = AffineForm(
        form.coeffs[: circuit.input_ct]
        + [0] * 3
        + form.coeffs[circuit.input_ct :]
        + [0] * circuit.mul_ct,
        form.constant,
    )
    assert len(newform.coeffs) == circuit.input_ct + 3 + 2 * circuit.mul_ct
    return newform


def calculate_fg_form(circuit, wire, challenge, gf):
    # Construct form for each mul-gate.
    forms = [construct_affine_form(g, circuit, wire) for g in circuit.mul_gates()]
    # Assert that form is consistent with length of z-vector.
    forms = [convert_to_ac20(f, circuit) for f in forms]

    lagr_range = range(circuit.mul_ct + 1)
    lagr_vect = lagrange(gf, lagr_range, challenge)

    # Construct form corresponding to poly evaluation at c, in coefficients of z-vector
    form = AffineForm([0] * circuit.input_ct + [0, 0, 0] + [0] * 2 * circuit.mul_ct, 0)
    form.coeffs[circuit.input_ct + wire] = 1 * lagr_vect[0]
    form += sum(forms[j] * l_j for j, l_j in enumerate(lagr_vect[1:]))
    return form


def calculate_h_form(circuit, challenge, gf):
    lagr_range = range(2 * circuit.mul_ct + 1)
    lagr_vect = lagrange(gf, lagr_range, challenge)
    form = LinearForm([0] * circuit.input_ct + [0] * 2 + lagr_vect)
    return form


def calculate_circuit_forms(circuit):
    forms = [
        construct_affine_form(circuit.gates[gate_ix], circuit, None)
        for gate_ix in circuit.output_gates
    ]
    return forms


def lagrange(gf, lagr_range, c):
    return verifiable_mpc.ac20.mpc_ac20._recombination_vectors(gf, lagr_range, (c,))[0]


# TODOs:
# 2. Add __sub__ gates (instead of -1*var + ...)
# 2. Name auxiliary inputs/witnesses? Witnesses and inputs are treated equally. However, witness can be created after inputs, e.g. for zero-test. Should we treat them separately?
# 3. circuit.output_gates list and circuit.out_gates() method: name consistently and unambiguously

# Notes
# https://stackoverflow.com/questions/2673651/inheritance-from-str-or-int
# https://stackoverflow.com/questions/20242479/printing-a-tree-data-structure-in-python
# https://github.com/egonSchiele/grokking_algorithms/blob/master/06_breadth-first_search/python/01_breadth-first_search.py
