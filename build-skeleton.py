#!/usr/bin/env python3

# %%
# imports
import networkx as nx
from enum import Enum, auto
import itertools as it
import random
from collections import defaultdict, Counter
from pprint import pprint
import matplotlib.pyplot as plt
from functools import reduce
import copy
import numpy as np
import math
import hashlib
import concurrent.futures

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
            self.c = []
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
        
    def wires(self):
        return [self.a, *sorted(self.c)]

    # For NOT and CNOT, default arguments give proper results
    def eval(self, a_val, c1_val=True, c2_val=True):
        return a_val ^ (c1_val & c2_val)

    def __str__(self, show_id=False):
        name = self.type.name if not show_id else f"{self.type.name}#{id(self)}"
        if self.type == GateT.NOT:
            return f"{self.type.name}[{self.a}]"
        return f"{name}[{self.a}; {' '.join(map(str, self.c))}]"
    
    def compact(self):
        short_type = {GateT.NOT: 'N', GateT.CNOT: 'C', GateT.CCNOT: 'T'}[self.type]
        return f"{short_type}[{self.a},{' '.join(map(str, self.c))}]"

    # Define a lexicographical ordering over Gates
    def __lt__(self, other, check_id=True):
        if len(self.c) < len(other.c):
            return True
        
        if len(self.c) > len(other.c):
            return False
        
        # same type. compare wires, and possibly IDs, if requested
        self_id = id(self) if check_id else None
        other_id = id(other) if check_id else None
        
        # use built-in tuple comparison
        return (self.a, *self.c, self_id) < (other.a, *other.c, other_id)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.type, self.a, tuple(self.c)))

    # def graph(self):
    #     """Generate the skeleton graph of this circuit."""


class CanonicalGate(Gate):
    def __init__(self, g):
        super().__init__(g.a, *sorted(g.c))

    def __eq__(self, other):
        return self.type == other.type and self.a == other.a and self.c == other.c

    def __lt__(self, other):
        r = super().__lt__(other, check_id=False)
        # print(f"{self} {'<' if r else '>='} {other}")
        return r

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
    gates = []
    def __init__(self, gates, wires=None):
        self.gates = []
        self.extend(gates)
        self.hash = hashlib.sha256(usedforsecurity=False)
        self.permdig = None

        if wires:
            self.wires = max(self.wires, wires)

    def extend(self, gates):
        self.gates.extend(copy.deepcopy(g) for g in gates)
        # zero wire is allowed
        self.wires = 1 + max(g.bottom() for g in gates)

    def __len__(self):
        return len(self.gates)
    
    def eval(self, input):
        assert(len(input) == self.wires)

        state = input[:]
        for g in self.gates:
            state[g.a] = g.eval(state[g.a], *(state[c] for c in g.c))

        return state

    def perm(self, digest=False):
        """return the permutation computed by this circuit, or a digest"""
        if digest and self.permdig:
            return self.permdig
        
        p = []
        for i in range(2 ** self.wires):
            b = format(i, "0{}b".format(self.wires))
            out_list = self.eval(list(map(int, b)))

            if digest:
                self.hash.update(bytes(out_list))
            else:
                out_n = reduce(lambda a, b: (a << 1) | b, out_list, 0)
                p.append(out_n)

        if digest:
            self.permdig = self.hash.hexdigest()
            return self.permdig
        else:
            return p
        
    def fingerprint(self):
        return self.perm(digest=True)

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
    
    def __repr__(self):
        return "{" + ' '.join(g.compact() for g in self.gates) + "}"

# %%
# A random circuit.


def random_circuit(n, m):
b = list(base_perms(n))
# deepcopy required here so each gate is "unique
return Circuit([random.choice(b) for _ in range(m)])


n = 8
m = 8
c = random_circuit(n, m)
print(c)

print(c.perm(digest=True))
print(c.perm(digest=False))

# %%
class PermutationDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._next_value = 0
    
    def __missing__(self, k):
        while self._next_value in self.values():
            self._next_value += 1
        self[k] = self._next_value
        return self._next_value

# %%
class SkeletonGraph(nx.DiGraph):
    """A lightweight class to represent skeleton graphs."""

    # edges do not need labels
    edge_dict = {}

    def from_circuit(ckt: Circuit):
        active_set = defaultdict(list)
        control_set = defaultdict(list)

        graph = {g: set() for g in ckt.gates}

        for g in ckt.gates:
            # get control pins on this active pin
            back = control_set[g.a]
            for h in filter(None, back):
                graph[h].update([g])
            
            control_set[g.a].append(None)

            if len(active_set[g.a]) > 0 and active_set[g.a][-1] is None:
                active_set[g.a] = [g]
            else:
                active_set[g.a].append(g)

            for c in g.c:
                # get active pins on each control pin
                for anode in filter(None, active_set.get(c, [])):
                    graph[anode].update([g])

                if len(control_set[c]) > 0 and control_set[c][-1] is None:
                    control_set[c] = [g]
                else:
                    control_set[c].append(g)
                active_set[c].append(None)

        # pprint(graph)
        return SkeletonGraph(graph)

    def single_edge_dict(self):
        return self.edge_dict

    edge_attr_dict_factory = single_edge_dict

    def canonical(self):
        topo = nx.topological_generations(self)

        wire_perm = PermutationDict()

        new_topo = [[]]

        for level in topo:
            ls = sorted(level)
            for g in ls:
                w = list(map(lambda p: wire_perm[p], g.wires()))
                # print(g.wires(), "->", w)
                new_topo[-1].append(CanonicalGate(Gate(*w)))
            new_topo.append([])

        # pprint(wire_perm)
        return tuple(tuple(lvl) for lvl in filter(None, new_topo))

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
g = SkeletonGraph.from_circuit(c)
print(g)

plt.figure(figsize=(6, 4))
g.draw()
plt.show()

for level in g.canonical():
    for cg in level:
        print(f"{str(cg):>15}", end=' ')
    print()


# %%
def all_circuits(n, m):
    b = base_perms(n)
    yield from map(lambda g: Circuit(g, wires=n), it.product(b, repeat=m))


# %% [md]
### Testing
# `4,4` takes ~ a minute and uses 2 GB RAM. Can be optimized.
n = 3
m = 3
K = all_circuits(n, m)
# %%

G = Counter(map(lambda c: SkeletonGraph.from_circuit(c).canonical(), K))

num_ckt = sum(G.values())

# %%
print(f"{len(G)} canonical; {num_ckt} circuits ({len(G) / num_ckt * 100:.1f}%)")

# %% [md]
#### Histogram
# pprint(G)

# Why is this not consistent?? Order of circuits shouldn't matter...
p = 0
# for k, v in sorted(cc.items(), key=lambda k: k[1], reverse=True):
#     print(f"{v:3} canonical circuits had {k:2} candidates each")
#     p += k * v

v1 = list(G.values())
# plt.figure(figsize=(8, 4))
# plt.hist(v1, bins=np.arange(1, max(v1) + 2), alpha=0.7)
# plt.show()
# v2 = list(G2.values())
# plt.hist(v2, bins=np.arange(1, max(v2) + 2), alpha=0.7)

# unique = len([k for k, c in G.items() if c == 1])
# print(f"{n=} {m=}: {unique} skeletons have a unique canonicalization")

G
# %%
def compact_repr(T):
    return ' '.join(''.join(g.compact() for g in L) for L in T)

for k, v in G.most_common():
    print(f"{compact_repr(k)}\t{v} circuits")

for m, ckt_c in Counter(G.values()).most_common():
    print(f"{ckt_c} skeleton graphs have circuit-multiplicity {m}")
# %%

GC = defaultdict(list)
count = Counter()
for c in K:
    canon = SkeletonGraph.from_circuit(c).canonical()
    GC[canon].append(c)
    count[canon] += 1
# %%

for k, v in count.most_common():
    print(k, v)
    for ckt in GC[k]:
        print(ckt)


# %%

c = Circuit([Gate(0), Gate(1, 0), Gate(1), Gate(2, 1)])
print(c)
# %%
g = SkeletonGraph.from_circuit(c)
# %%
g.draw()
# %%


len(G) / num_ckt
# %% [md]
### Preprocess all circuits
#
# For each circuit sampled, build the canonicalized skeleton graph and determine
# the permutation it computes. Bucket skeleton graphs by permutation.

class Permutation(tuple):
    def __init__(self, t):
        self.t = t

    def __str__(self):
        # TODO: cycle notation
        return "[" + " ".join(
            map(
                lambda x: '.' if x[0] == x[1] else str(x[1]),
                zip(range(len(self.t)), self.t)
            )) + "]"

class SkeletonCache():
    def __init__(self, n, m):
        self.n = n
        self.m = m
        print(f"Building Skeleton Cache with {n=}, {m=}")

        self.b = list(base_perms(n))
        # TODO: up to reordering
        self.n_perms = math.factorial(2 ** n) // 2
        print(f"{self.n_perms} perms on {n} wires. b={len(self.b)}; {len(self.b) ** m} circuits.")

    def build(self):
        K = all_circuits(self.n, self.m)
        G = defaultdict(lambda: defaultdict(list))

        canon2perm = {}

        def ckt_work(ckt):
            return Permutation(ckt.perm()), SkeletonGraph.from_circuit(ckt).canonical()

        # with concurrent.futures.ProcessPoolExecutor() as exec:
        #     for ckt, (perm, graph) in zip(K, exec.map(ckt_work, K)):
        #         G[perm][graph].append(ckt)
        

        for ckt in K:
            perm, canon = ckt_work(ckt)
            # get the canonical perm
            cperm = canon2perm.setdefault(canon, perm)

            G[cperm][canon].append(ckt)
        
        # graphs = sum(len(v) for _, v in G.items())
        # print(f"{graphs} graphs generated {len(G)} perms ({pct*100:.1f}%)")
        self.G = G
        return G
    
    def stats(self):
        num_graphs = sum(len(v) for v in self.G.values())
        num_perms = len(G)

        iden_perm = Permutation(range(2 ** self.n))
        iden_graphs = G[iden_perm]
        num_id_graphs = len(iden_graphs)
        num_id_ckts = sum(len(v) for v in iden_graphs.values())

        uniq_graph = sum(1 if len(v) == 1 else 0 for v in G.values())
        uniq_ckt = sum(1 if (sum(len(c) for c in v) == 1) else 0 for v in G.values())
    
        print(f"{len(self.b) ** self.m} circuits => {num_graphs} uniq graphs => {num_perms} uniq perms")
        print(f"  {num_id_ckts} id ckts => {num_id_graphs} id graphs")
        print(f"  {uniq_graph} perms w/ a uniq graph")
        print(f"  {uniq_ckt} perms w/ a uniq circuit")

    def print(self):
        for p, gs in self.G.items():
            title = f"Perm {p}: {len(gs)} graphs, {sum(len(k) for k in gs.values())} circuits"
            print("\n" + "=" * len(title))
            print(title)

            for i, (g, ckt) in enumerate(gs.items()):
                print(f"#{1+i} Graph {g}: {len(ckt)} circuits.")

                for k in ckt:
                    print(k)
                    break


# %%
cache = SkeletonCache(n=4, m=3)
G = cache.build()

cache.stats()

# %%
cache.print()

# %%
pprint(G)
# %%
count_ckt = [sum(v.values()) for _, v in G.items()]
count_graph = [len(v) for _, v in G.items()]
# %%
plt.hist(count_ckt, bins=np.arange(1, max(count_ckt) + 2), alpha=0.5)
ax = plt.twinx()
ax.hist(count_graph, bins=np.arange(1, max(count_graph) + 2), color='orange', alpha=0.5)
plt.xlim(0, 50)

# %%
s = sorted(G, key=lambda k:sum(G[k].values()), reverse=True)
for k in s:
    print(k, sum(G[k].values()))
# %%
pprint(G)
# %%


sum(1 if len(v) == 1 else 0 for v in G.values())
# %%
len(G)
# %%


def second(x):
    return tuple(map(lambda y: y[1], x))

G = map(SkeletonGraph.from_circuit, all_circuits(3, 3))

dv = map(lambda g: (tuple(zip(second(g.in_degree), second(g.out_degree)))), G)
s = set(dv)
print(f"{len(s)} unique degree assignments")
pprint(s)

# %%
c = random_circuit(8, 8)
print(c)

# %%
from circuit_lib import SkeletonGraph as sg
c = random_circuit(8, 8)
g = sg.from_circuit(c, True)


## TODO FROM RENE: color edges

# %%
print(c)
g.draw()
# %%

# %%
