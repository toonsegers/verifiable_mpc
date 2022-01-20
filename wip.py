"""WIP: Plot circuit graphs with networkx. *unfinished*
"""
import pprint
from mpyc.finfields import GF
import verifiable_mpc.ac20.circuit_sat_cb as cs
import verifiable_mpc.ac20.circuit_builder as cb
from mpyc.fingroups import QuadraticResidues

import matplotlib.pyplot as plt
import networkx as nx


pp = pprint.PrettyPrinter(indent=4)


def save_circuit(circuit):
    # Maybe: Use list returned by this method for print_circuit (DRY)
    ret = []
    for gate in circuit.out_gates():
        ret += save_out_gate(circuit, gate)
    return ret


def save_out_gate(circuit, gate, level=0):
    ret = [{
        "level": level,
        "in": [i.name if isinstance(i, cb.CircuitVar) else i for i in gate.inputs], 
        "out": gate.output.name if isinstance(gate.output, cb.CircuitVar) else gate.output,
        "op": gate.op,
        "index": gate.index
        }]
    for child in circuit.children(gate):
        ret += save_out_gate(circuit, child, level + 1)
    return ret


def main(pivot_choice, n):
    print("Pivot selected: ", pivot_choice)

    group = QuadraticResidues(l=1024)
    gf = GF(modulus=group.order)

    circuit = cb.Circuit()
    # Using field elements
    # b = cb.CircuitVar(gf(1), circuit, "b")
    # c = cb.CircuitVar(gf(2), circuit, "c")
    # Using integers
    b = cb.CircuitVar(1, circuit, "b")
    c = cb.CircuitVar(2, circuit, "c")

    # Example: Mini circuit
    # d = c + c 
    # e = d * c  
    # e.label_output("e")

    # Example: Not-equal circuit
    g = b != c
    g.label_output("g")

    # Example: larger circuit with >n mul-gates
    # d = c + c + c * c + c * c * 1 + 1 + b 
    # e = d*d + c**n + 10
    # f = d*c + e
    # f.label_output("f")
    # g = f != 100
    # g.label_output("g")
    # h = g >= 10  # Note: comparison only works for integers
    # h.label_output("h")

    print(circuit)

    circ_lst = save_circuit(circuit)
    print(circ_lst)

    # assumes binary tree (fan-in 2) just as AC20
    notflat = [[
            ((i["in"][0], i["level"]), (i["op"], i["index"], i["level"])), 
            ((i["in"][1], i["level"]), (i["op"], i["index"], i["level"])), 
            ((i["op"], i["index"], i["level"]), (i["out"], i["level"]-1)), 
            ] for i in circ_lst]

    edges = [j for sub in notflat for j in sub]
    print(edges)

    G = nx.MultiDiGraph()
    G.add_edges_from(edges)
    print("Edges:")
    print([e for e in G.edges])
    print("Nodes:")
    print([n for n in G.nodes])

    print("Is tree/arborescence/forest:")
    print(nx.is_tree(G))
    print(nx.is_arborescence(G))
    print(nx.is_forest(G))

    pos = nx.nx_agraph.graphviz_layout(G, prog="dot", args="")
    plt.figure(figsize=(20, 10))
    nx.draw(G, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=True)
    plt.axis("equal")
    plt.savefig('circuit.png')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="roughly number of multiplications")
    parser.set_defaults(n=3)
    args = parser.parse_args()

    verification = main(cs.PivotChoice.compressed, args.n)



# Graveyard
    # nx.draw(G, with_labels=True)
    # plt.savefig('hierarchy.png')

    # # assumes binary tree (fan-in 2) just as AC20

    # # notflat = [[(t["in"][0], t["out"]), (t["in"][1], t["out"])] for i, t in enumerate(circ_lst)]
    # # notflat = [[(i["in"][0], i["out"]), (i["in"][1], i["out"])] for i in circ_lst]
    # # edges = [j for sub in [[(i["in"][0], i["out"]), (i["in"][1], i["out"])] for i in circ_lst] for j in sub]

    # # TODO: make gates part of the graph
    # notflat = [[
    #         ((i["in"][0],0,i["level"]), (i["out"],0,i["level"]-1)), 
    #         ((i["in"][1],1,i["level"]), (i["out"],0,i["level"]-1)), 
    #         ] for i in circ_lst]
    # edges = [j for sub in notflat for j in sub]
    # print(edges)

    # G = nx.MultiDiGraph()
    # G.add_edges_from(edges)
    # print("Edges:")
    # print([e for e in G.edges])
    # print("Nodes:")
    # print([n for n in G.nodes])

    # print("Is tree/arborescence/forest:")
    # print(nx.is_tree(G))
    # print(nx.is_arborescence(G))
    # print(nx.is_forest(G))

    # nx.draw(G, with_labels=True)
    # plt.savefig('hierarchy.png')

    # pos = nx.nx_agraph.graphviz_layout(G, prog="dot", args="")
    # plt.figure(figsize=(100, 100))
    # nx.draw(G, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=True)
    # plt.axis("equal")
    # plt.savefig('hierarchy2.png')