"""Demonstrate circuit_builder module.

Initialize a circuit using the Circuit class. 
The Circuit class tracks gates and variables corresponding to a circuit.
A circuit is constructed by using instances of CircuitVar in arithmetic.
"""
import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
    
from mpyc.finfields import GF 
from mpyc.runtime import mpc
import verifiable_mpc.ac20.circuit_builder as cb


if __name__ == "__main__":
    circuit = cb.Circuit()
    
    # Initialize field with modulus 2^31-1; https://en.wikipedia.org/wiki/2,147,483,647
    gf = GF(2_147_483_647)  

    # Initialize input variables b, c.
    b = cb.CircuitVar(gf(2), circuit, "b")
    c = cb.CircuitVar(gf(2), circuit, "c")

    # Construct arithmetic circuit using pure Python.
    # d = 10 - c 
    # d = c - 10
    d = c + c + c * c + c * c * 1 + 1 + b 
    e = d * d + c + 10
    f = d * c + e
    # f = b * b + c - c
    # f = b * b 

    # Label the variables that represent output gates.    
    f.label_output("f")
    g = f + 100
    g.label_output("g")
    print(f"Output gates: {f=} and {g=}")

    # Print attributes of the circuit up to this point (circuits can extended).
    print("Circuit attributes:")
    print(f"{circuit.mul_ct=}")
    print(f"{circuit.add_ct=}")
    print(f"{circuit.input_ct=}")
    print(f"{circuit.output_ct=}")
    print("String representation of circuit:")
    print(cb.print_circuit(circuit))
    print("CircuitVars: ", circuit.circuitvars)
    print("Input gates (indexes): ", circuit.input_gates)

    # Evaluate the circuit for given inputs:
    inputs = [gf(2), gf(2)]
    print(f"Evaluate circuit for {inputs=}: {circuit(inputs)}")

    # print(f"Evaluate circuit up to and including a specified gate:")
    # for g in circuit.gates:
    #     y = circuit.eval(inputs, g)
    #     print(f"{y=} for {str(g)}")

    print("Construct objects for AC20 proof system.")
    alpha, beta, gamma = circuit.multiplication_triples(inputs)
    f_poly, g_poly, h_poly = cb.calculate_fgh_polys(alpha, beta, None, gf) 
    challenge = 15
    f_form = cb.calculate_fg_form(circuit, wire=0, challenge=challenge, gf=gf)
    g_form = cb.calculate_fg_form(circuit, wire=1, challenge=challenge, gf=gf)
    h_form = cb.calculate_h_form(circuit, challenge, gf)
    print(f"{f_form=}")
    print(f"{g_form=}")
    print(f"{h_form=}")

    # TODO: construct z_vector without construction of f_, g_, h_poly (see mpc_ac20.py)
    z_vector = inputs + [f_poly(0), g_poly(0)] + [h_poly(i) for i in range(2 * circuit.mul_ct + 1)]

    assert f_form(z_vector) * g_form(z_vector) == h_form(z_vector)
    assert f_poly(challenge) * g_poly(challenge) == h_poly(challenge)
    assert h_form(z_vector) == h_poly(challenge)

    print("Construct form(s) corresponding to circuit output(s), and evaluate up to and including that gate.")
    # See also Section 5.2 and 6.1 of 'Pinocchio-Based Adaptive zk-SNARKs and Secure/Correct Adaptive Function Evaluation' [Vee17]
    # Link to [Vee17]: https://eprint.iacr.org/2017/013.pdf
    forms = cb.calculate_circuit_forms(circuit)
    print("Circuit output forms: ", forms)
    print([form(inputs+gamma) for form in forms])


    print("*** Demo gadgets ***")
    print("Demonstrate != 0 gadget, particularly b = (a != 0) ? 1 : 0")
    # Start with new circuit
    circuit = cb.Circuit()
    # Initialize input variables b, c.
    a = cb.CircuitVar(gf(3), circuit, "a") # Input variable, in this case != 0
    print(f"Input a = {a}")
    b = a != 0  # b = (a != 0) ? 1 : 0
    b.label_output("b")
    print("b=", b)
    # Circuit builder writes equations a · c = b, a · (1 − b) = 0  
    print(circuit)
    a = gf(3)
    # Retrieve inputs (including auxiliary inputs) that were recorded during circuit construction.
    inputs = circuit.initial_inputs()
    print(f"Evaluate circuit for {inputs=}: {circuit(inputs)}")
    
    # Demo a <= b gadget, start with new circuit.
    print("Demonstrate <= gadget.")
    circuit = cb.Circuit()
    secint = mpc.SecInt()
    a = cb.CircuitVar(secint(5), circuit, "b")
    b = cb.CircuitVar(secint(3), circuit, "c")
    c = a <= b
    c.label_output("c")
    # print(cb.print_out_gate(circuit, circuit.out_gates()[0]))
    # print(circuit)



