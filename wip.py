"""WIP: Plot circuit graphs with networkx. *unfinished*
"""

import random
import pprint
import sys, os

project_root = sys.path.append(os.path.abspath(".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from mpyc.finfields import GF

import verifiable_mpc.ac20.circuit_sat_cb as cs
import verifiable_mpc.ac20.circuit_builder as cb
from sec_groups.fingroups import QuadraticResidue, EllipticCurve
import sec_groups.ellcurves as ell
from sec_groups.tools.find_primes import find_safe_primes

import matplotlib.pyplot as plt
import networkx as nx


pp = pprint.PrettyPrinter(indent=4)

PIVOT = cs.PivotChoice.compressed  
GROUP = "QR"




def hierarchy_pos(G, root=None, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5):

    '''
    From Joel's answer at https://stackoverflow.com/a/29597209/2966723.  
    Licensed under Creative Commons Attribution-Share Alike 
    
    If the graph is a tree this will return the positions to plot this in a 
    hierarchical layout.
    
    G: the graph (must be a tree)
    
    root: the root node of current branch 
    - if the tree is directed and this is not given, 
      the root will be found and used
    - if the tree is directed and this is given, then 
      the positions will be just for the descendants of this node.
    - if the tree is undirected and not given, 
      then a random choice will be used.
    
    width: horizontal space allocated for this branch - avoids overlap with other branches
    
    vert_gap: gap between levels of hierarchy
    
    vert_loc: vertical location of root
    
    xcenter: horizontal location of root
    '''
    if not nx.is_tree(G):
        raise TypeError('cannot use hierarchy_pos on a graph that is not a tree')

    if root is None:
        if isinstance(G, nx.DiGraph):
            root = next(iter(nx.topological_sort(G)))  #allows back compatibility with nx version 1.11
        else:
            root = random.choice(list(G.nodes))

    def _hierarchy_pos(G, root, width=1., vert_gap = 0.2, vert_loc = 0, xcenter = 0.5, pos = None, parent = None):
        '''
        see hierarchy_pos docstring for most arguments

        pos: a dict saying where all nodes go if they have been assigned
        parent: parent of this branch. - only affects it if non-directed

        '''
    
        if pos is None:
            pos = {root:(xcenter,vert_loc)}
        else:
            pos[root] = (xcenter, vert_loc)
        children = list(G.neighbors(root))
        if not isinstance(G, nx.DiGraph) and parent is not None:
            children.remove(parent)  
        if len(children)!=0:
            dx = width/len(children) 
            nextx = xcenter - width/2 - dx/2
            for child in children:
                nextx += dx
                pos = _hierarchy_pos(G,child, width = dx, vert_gap = vert_gap, 
                                    vert_loc = vert_loc-vert_gap, xcenter=nextx,
                                    pos=pos, parent = root)
        return pos

            
    return _hierarchy_pos(G, root, width, vert_gap, vert_loc, xcenter)


def main(pivot_choice, group_choice, n):
    print("Pivot selected: ", pivot_choice)

    if pivot_choice == cs.PivotChoice.koe:
        # TODO: improve syntax for passing two groups to create_generators
        group1 = EllipticCurve(ell.BN256, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        group2 = EllipticCurve(ell.BN256_TWIST, ell.WEI_HOM_PROJ, ell.Weierstr_HomProj_Arithm)
        group1.is_additive = False
        group1.is_multiplicative = True
        group2.is_additive = False
        group2.is_multiplicative = True
        group = [group1, group2]
        order = group1.order
        gf = GF(modulus=order)
    elif GROUP == "Elliptic":
        group = EllipticCurve(ell.ED25519, ell.ED_HOM_PROJ, ell.Edwards_HomProj_Arithm)
        group.is_additive = False
        group.is_multiplicative = True
        gf = GF(modulus=group.order)
    elif GROUP == "QR":
        order, modulus = find_safe_primes(1024)
        group = QuadraticResidue(modulus=modulus)
        gf = GF(modulus=group.order)

    circuit = cb.Circuit()
    # b = cb.CircuitVar(gf(1), circuit, "b")
    # c = cb.CircuitVar(gf(2), circuit, "c")
    b = cb.CircuitVar(1, circuit, "b")
    c = cb.CircuitVar(2, circuit, "c")

    d = c + c + c * c + c * c * 1 + 1 + b 
    e = d*d + c**n + 10
    f = d*c + e
    f.label_output("f")
    g = f != 100
    g.label_output("g")
    h = g >= 10  # Note: comparison only works for integers
    h.label_output("h")


    def collect_edges(circuit, gate, level=0):
        edges = []
        nodes = [(gate.index, {"op": gate.op})]
        for child in circuit.children(gate):
            edges += [(child.index, gate.index)]
            e, n = collect_edges(circuit, child, level + 1)
            edges += e
            nodes += n
        return edges, nodes

    edges = []
    nodes = []
    for gate in circuit.out_gates():
        e, n = collect_edges(circuit, gate)
        edges += e 
        nodes += n
    # edges, nodes = collect_edges(circuit, circuit.out_gates()[-1])

    G=nx.Graph()
    G.add_edges_from(edges)
    print([e for e in G.edges])
    print([n for n in G.nodes])
    # pos = hierarchy_pos(G,1)    
    # pos = nx.nx_agraph.graphviz_layout(G, prog="twopi", args="")
    # plt.figure(figsize=(8, 8))
    # nx.draw(G, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=False)
    # # nx.draw(G, pos=pos, with_labels=True)
    # nx.draw(G, with_labels=True, font_weight='bold')
    plt.savefig('hierarchy.png')




    # x = circuit.initial_inputs()
    # # Check if resulting commitment vector is of appropriate length.
    # check, padding, g_length = cs.check_input_length_power_of_2(x, circuit)
    # # Add unused variables to pad the length of the commitment vector to power of 2 minus 1.
    # unused = [cb.CircuitVar(0, circuit, "unused_"+str(i)) for i in range(padding)]
    # x = circuit.initial_inputs()
    # print("Length of input vector including auxiliary inputs (witnesses for special gates): ", len(x))
    # print("Length of commitment vector: ", g_length)

    # generators = cs.create_generators(g_length, pivot_choice, group, progress_bar=True)
    # print("Generators created/trusted setup done.")

    # print("Start non-interactive circuit satisfiability proof with compressed pivot. ")
    # proof = cs.circuit_sat_prover(generators, circuit, x, gf, pivot_choice)
    # # print("Proof:")
    # # pp.pprint(proof)
    # print("Start verification.")
    # verification = cs.circuit_sat_verifier(proof, generators, circuit, gf, pivot_choice)
    # print("Verification checks: ")
    # pp.pprint(verification)

    # return verification



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", type=int, help="roughly number of multiplications")
    parser.add_argument("--elliptic", action="store_true", 
                        help="use elliptic curve groups (default QR groups)")
    parser.add_argument('--basic', action='store_true',
                        help='use basic pivot (not the compressed pivot)')
    parser.add_argument('--koe', action='store_true',
                        help='use pivot based on Knowledge-of-Exponent assumption and BN256 curves')
    parser.set_defaults(n=3)
    args = parser.parse_args()
    if args.elliptic:
        GROUP = "Elliptic"
    elif args.basic:
        PIVOT = cs.PivotChoice.pivot
    elif args.koe:
        PIVOT = cs.PivotChoice.koe

    verification = main(PIVOT, GROUP, args.n)



# TODO:
# TS-1: Make interface of create_generators(,,group,) more intuitive when working with KOE
# TS-1: Remove L from proof; non-interactive nullity proof does not have to pass it to verifier.
