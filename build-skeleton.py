#!/usr/bin/env python3

# %%
# imports
import networkx as nx
from enum import Enum, auto
import itertools as it
import random
from collections import defaultdict
from pprint import pprint
import matplotlib.pyplot as plt
import functools
import copy

# %% [md]
# # Skeleton Graphs
#
# This python file (/ interactive notebook) will construct all unique skeleton
# graphs for reversible circuits of a given size (wires x gates). For our base
# gates we are using the (2-input) Toffoli gate, (1-input) CNOT gate, and
# NOT gate.
#
# For a reversible circuit with $n$ wires and $m$ gates, each of which may be
# one of $b=o(n^3)$ base permutations, there are $b^m$ possible circuits.
# However, many of these will _not_ have unique skeleton graphs. Analytically
# computing the number of such unique graphs is nontrivial because it would
# require determining the automorphisms over the permutation group for each
# graph and invoking Burnside's Lemma.
#
# Skeleton graphs do not take into account permutations of the wires.
# %%
# Compute base permutations


class GateT(Enum):
    NOT = auto()
    CNOT = auto()
    CCNOT = auto()


class Gate:
    """A general CxNOT gate."""

    c = []

    def __init__(self, a, c1=None, c2=None):
        if c1 is None:
            # NOT gate
            self.a = a
            self.type = GateT.NOT
        elif c2 is None:
            # CNOT gate
            self.a = a
            self.c = [c1]
            self.type = GateT.CNOT
        else:
            # CCNOT gate
            self.a = a
            self.c = [c1, c2]
            self.type = GateT.CCNOT

    def top(self):
        if self.type == GateT.NOT:
            return self.a
        else:
            return min(self.a, min(self.c))

    def bottom(self):
        if self.type == GateT.NOT:
            return self.a
        else:
            return max(self.a, max(self.c))

    # For NOT and CNOT, default arguments give proper results
    def eval(self, a_val, c1_val=True, c2_val=True):
        return a_val ^ (c1_val & c2_val)

    def __str__(self, show_id=False):
        name = self.type.name if not show_id else f"{self.type.name}#{id(self)}"
        if self.type == GateT.NOT:
            return f"{self.type.name}[{self.a}]"
        return f"{name}[{self.a}; {' '.join(map(str, self.c))}]"

    # Define a lexicographical ordering over Gates
    def __lt__(self, other, check_id=True):
        if self.a < other.a:
            return True

        if len(self.c) > 0 and len(other.c) > 0:
            if self.c[0] < other.c[0]:
                return True
            if len(self.c) == 2 and len(other.c) == 2:
                if self.c[1] < other.c[1]:
                    return True

        if check_id:
            return id(self) < id(other)
        else:
            return False

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.type, self.a, tuple(self.c)))

    # def graph(self):
    #     """Generate the skeleton graph of this circuit."""


class CanonicalGate(Gate):
    def __init__(self, g):
        super().__init__(g.a, *g.c)

    def __eq__(self, other):
        return self.type == other.type and self.a == other.a and self.c == other.c

    def __hash__(self):
        return super().__hash__()


def base_perms(n):
    """Generates all base permutations over $n$ wires."""

    # Tof gates
    for active in range(n):
        other = (i for i in range(n) if i != active)
        yield from (Gate(active, c[0], c[1]) for c in it.combinations(other, 2))

    # CNOT gates
    for active in range(n):
        other = (i for i in range(n) if i != active)
        yield from (Gate(active, c) for c in other)

    # NOT gates
    yield from (Gate(active) for active in range(n))


# %% [md]
# Base perms over $n=4$ wires:
n = 4
b = list(base_perms(n))
print(f"{n=}: {len(b)} base perms")
for g in b:
    print(g)

# %%

# Box drawing characters
CKT_WIRE_CHAR = "─"
CKT_CONTROL_CHAR = "│"
CKT_DOT_CHAR = "●"


class Circuit:
    def __init__(self, gates):
        self.gates = []
        self.extend(gates)

    def extend(self, gates):
        self.gates.extend(copy.deepcopy(g) for g in gates)
        # zero wire is allowed
        self.wires = 1 + max(g.bottom() for g in gates)

    def __len__(self):
        return len(self.gates)

    def __str__(self):
        s = ""
        for w in range(self.wires):
            s += f"{w:>3} {CKT_WIRE_CHAR}"
            for g in self.gates:
                if w == g.a:
                    if g.type == GateT.NOT:
                        s += "[!]"
                    else:
                        s += "[&]"
                elif w in g.c:
                    s += CKT_WIRE_CHAR + CKT_DOT_CHAR + CKT_WIRE_CHAR
                elif w > g.top() and w < g.bottom():
                    s += CKT_WIRE_CHAR + CKT_CONTROL_CHAR + CKT_WIRE_CHAR
                else:
                    s += 3 * CKT_WIRE_CHAR

                s += CKT_WIRE_CHAR

            s += "\n"

            s += " " * (4 if w == self.wires - 1 else 5)

            for m, g in enumerate(self.gates):
                if w == self.wires - 1:
                    s += f"{m:3} "
                elif w >= g.top() and w < g.bottom():
                    s += " " + CKT_CONTROL_CHAR + "  "
                else:
                    s += " " * 4

            s += "\n"

        return s


# %%
# A random circuit.


def random_circuit(n, m):
    b = list(base_perms(n))
    # deepcopy required here so each gate is "unique
    return Circuit([random.choice(b) for _ in range(m)])


n = 5
m = 8
c = random_circuit(n, m)
print(c)


# %%
class SkeletonGraph(nx.DiGraph):
    """A lightweight class to represent skeleton graphs."""

    # edges do not need labels
    edge_dict = {}

    def from_circuit(ckt: Circuit):
        active_set = {}
        control_set = defaultdict(list)

        graph = {g: set() for g in ckt.gates}

        for g in ckt.gates:
            # get control pins on this active pin
            back = control_set[g.a]
            for h in back:
                graph[h].update([g])
            active_set[g.a] = g
            control_set[g.a] = []

            for c in g.c:
                # get active pins on each control pin
                anode = active_set.get(c, None)
                if anode is not None:
                    graph[anode].update([g])

                control_set[c].append(g)

        # pprint(graph)
        return SkeletonGraph(graph)

    def single_edge_dict(self):
        return self.edge_dict

    edge_attr_dict_factory = single_edge_dict

    def canonical(self):
        topo = nx.topological_generations(self)

        return tuple(tuple(CanonicalGate(g) for g in sorted(level)) for level in topo)

    def draw(self):
        gen = list(nx.topological_generations(self))

        for layer, nodes in enumerate(gen):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                g.nodes[node]["layer"] = layer

        pos = nx.multipartite_layout(g, subset_key="layer")

        nx.draw(
            g,
            pos=pos,
            with_labels=False,
            connectionstyle="arc3,rad=0.4",
            node_color="black",
            edge_color=[random.choice(plt.cm.Dark2.colors) for _ in g.nodes],
            node_size=1200,
        )

        labels = {node: str(node).replace("[", "\n")[:-1] for node in g}
        nx.draw_networkx_labels(
            g, pos=pos, labels=labels, font_size=8, font_color="white", font_weight="bold"
        )

# %%
print(c)
g = SkeletonGraph.from_circuit(c)
print(g)

plt.figure(figsize=(6, 4))
g.draw()

for level in g.canonical():
    for cg in level:
        print(f"{str(cg):>15}", end=' ')
    print()


# %%
def all_circuits(n, m):
    b = base_perms(n)
    yield from map(Circuit, it.product(b, repeat=m))


# %%
K = list(all_circuits(3, 4))
# %%
G = set(map(lambda c: SkeletonGraph.from_circuit(c).canonical(), K))
# %%

print(len(G), len(K), f"{len(G) / len(K) * 100:.1f}%")

# %%
