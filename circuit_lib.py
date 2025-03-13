from collections import defaultdict, Counter
from enum import Enum, auto
from functools import reduce, partial
import copy
import numpy as np
import itertools as it
import networkx as nx
import random
import matplotlib.pyplot as plt
import math
from rich.progress import (
    Progress,
    BarColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
    MofNCompleteColumn,
    TransferSpeedColumn,
)
from rich.pretty import pprint
from rich import print
import multiprocessing as mp
import os

class GateT(Enum):
    # A regular CNOT gate.
    # a' = a ^ b
    CNOT   = "C"
    # An inverted CNOT gate
    # a' = ! (a ^ b)
    nCNOT  = "nC"
    # Toffoli gate
    # a' = a ^ (b & c)
    TOF  = "T"
    # Inverted Toffoli gate
    # a' = a  ^ !(b & c)
    nTOF = "nT"


class Gate:
    """A general CxNOT gate."""

    c = []
    inverted = False

    def __init__(self, a, c1, c2=None, inverted=False):
        assert type(a) is int
        self.inverted = inverted
        
        if c2 is None:
            # CNOT gate
            assert type(c1) is int
            self.a = a
            self.c = (c1,)
            self.type = GateT.nCNOT if inverted else GateT.CNOT
        else:
            # TOF gate
            assert type(c1) is int and type(c2) is int
            self.a = a
            self.c = (c1, c2)
            self.type = GateT.nTOF if inverted else GateT.TOF

    # extents of the gate
    def top(self):
        return min(self.a, min(self.c))

    def bottom(self):
        return max(self.a, max(self.c))

    def wires(self):
        return [self.a, *sorted(self.c)]

    # For CNOT, default argument gives proper results
    def eval(self, a_val, c1_val, c2_val=True):
        return a_val ^ self.inverted ^ (c1_val & c2_val)

    def __str__(self, show_id=False):
        name = self.type.value if not show_id else f"{self.type.value}#{id(self)}"
        return f"{name}[{self.a}; {' '.join(map(str, self.c))}]"

    def compact(self):
        short_type = self.type.value
        return f"{short_type}[{self.a},{' '.join(map(str, self.c))}]"

    # Define a lexicographical ordering over Gates
    def __lt__(self, other, check_id=True):
        # CNOTs before TOFs
        if len(self.c) < len(other.c):
            return True

        if len(self.c) > len(other.c):
            return False

        # same type. compare wires, and possibly IDs, if requested
        self_id = id(self) if check_id else None
        other_id = id(other) if check_id else None

        # use built-in tuple comparison
        # ID is used to disambiguate identically-connected gates
        return (self.inverted, self.a, *self.c, self_id) < (other.inverted, other.a, *other.c, other_id)

    def __repr__(self):
        return self.__str__()

    def __hash__(self):
        return hash((self.type, self.a, tuple(self.c)))
    
    def copy(self):
        if len(self.c) == 1:
            return Gate(self.a, self.c[0], None, self.inverted)
        else:
            return Gate(self.a, *self.c, self.inverted)


class CanonicalGate(Gate):
    def __eq__(self, other):
        return self.type == other.type and self.a == other.a and self.c == other.c

    def __lt__(self, other):
        r = super().__lt__(other, check_id=False)
        return r

    def __hash__(self):
        return super().__hash__()

def base_perms(n):
    """Generates all base permutations over $n$ wires."""

    # Tof gates
    for active in range(n):
        o1, o2 = it.tee((i for i in range(n) if i != active))
        # TOF
        yield from (Gate(active, c[0], c[1]) for c in it.combinations(o1, 2))
        # nTOF
        yield from (Gate(active, c[0], c[1], True) for c in it.combinations(o2, 2))

    # CNOT gates
    for active in range(n):
        o1, o2 = it.tee(i for i in range(n) if i != active)
        # CNOT
        yield from (Gate(active, c) for c in o1)
        # nCNOT
        yield from (Gate(active, c, None, True) for c in o2)

def zero_active_perms(n):
    b = base_perms(n)
    yield from filter(lambda g: g.a == 0, b)

class Circuit:
    gates = []

    def __init__(self, gates, wires=None):
        self.gates = []
        self.extend(gates)
        # self.hash = hashlib.sha256(usedforsecurity=False)
        # self.permdig = None

        if wires:
            self.wires = max(self.wires, wires)

    def extend(self, gates):
        self.gates.extend(g.copy() for g in gates)
        # zero wire is allowed
        self.wires = 1 + max(g.bottom() for g in gates)

    def __len__(self):
        return len(self.gates)

    def eval(self, input):
        assert len(input) == self.wires

        state = input[:]
        for g in self.gates:
            state[g.a] = g.eval(state[g.a], *(state[c] for c in g.c))

        return state

    def perm(self, digest=False):
        """return the permutation computed by this circuit, or a digest"""
        if digest and self.permdig:
            return self.permdig

        p = []
        for i in range(2**self.wires):
            b = format(i, "0{}b".format(self.wires))
            out_list = self.eval(list(map(int, b)))

            if digest:
                pass
                # self.hash.update(bytes(out_list))
            else:
                out_n = reduce(lambda a, b: (a << 1) | b, out_list, 0)
                p.append(out_n)

        if digest:
            # self.permdig = self.hash.hexdigest()
            return self.permdig
        else:
            return tuple(p)
        
    

    # def __eq__(self, other):
    #     if len(self) != len(other):
    #         return False

    #     return all(CanonicalGate(g1) == CanonicalGate(g2) for g1, g2 in zip(self.gates, other.gates))

    # def __hash__(self):
    #     return hash(tuple(CanonicalGate(g) for g in self.gates))

    def fingerprint(self):
        return self.perm(digest=True)

    def __str__(self):
        CKT_WIRE_CHAR = "─"
        CKT_CONTROL_CHAR = "│"
        CKT_DOT_CHAR = "●"
        s = ""
        for w in range(self.wires):
            s += f"{w:>3} {CKT_WIRE_CHAR}"
            for g in self.gates:
                if w == g.a:
                    if g.inverted:
                        s += "(!)"
                    else:
                        s += "( )"
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
        return "{" + " ".join(g.compact() for g in self.gates) + "}"


def all_circuits(n, m, zero_start=False):
    b = list(base_perms(n))

    if zero_start:
        z = list(zero_active_perms(n))
        itl = [z] + [b] * (m - 1)
        yield from it.product(*itl)
    else:
        yield from it.product(b, repeat=m)


def random_circuit(n, m):
    b = list(base_perms(n))
    return Circuit([random.choice(b) for _ in range(m)])

def all_unique(x):
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)

class PermutationDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._next_value = 0

    def __missing__(self, k):
        while self._next_value in self.values():
            self._next_value += 1
        self[k] = self._next_value
        return self._next_value


class NodeType(Enum):
    Wire = 'W'
    Gate = 'G'

    def __lt__(self, other):
        return self.value < other.value
    
    def __repr__(self):
        return self.value
    
class EdgeWeight(complex):
    def __lt__(self, y):
        s_re = self.real
        s_im = self.imag
        y_re = y.real
        y_im = y.imag
        
        s_norm = (abs(s_re) + abs(s_im), s_re, s_im)
        y_norm = (abs(y_re) + abs(y_im), y_re, y_im)
        return s_norm < y_norm

class SkeletonGraph(nx.DiGraph):
    """A lightweight class to represent skeleton graphs."""

    def from_circuit(ckt: Circuit, augmented=False):
        w = ckt.wires
        active_set = [[] for _ in range(w)]
        control_set = [[] for _ in range(w)]

        sk = SkeletonGraph()

        if augmented:
            sk.add_nodes_from(range(ckt.wires))

        sk.add_nodes_from(ckt.gates)

        for g in ckt.gates:
            # get control pins on this active pin
            back = control_set[g.a]
            for h in filter(None, back):
                # this is a control -> active collision
                sk.add_edge(h, g, wt=1j)

            if augmented:
                # also add the wire edges
                sk.add_edge(g, g.a, wt=1)
                for c_ in g.c:
                    sk.add_edge(g, c_, wt=1j)

            control_set[g.a].append(None)

            if len(active_set[g.a]) > 0 and active_set[g.a][-1] is None:
                active_set[g.a] = [g]
            else:
                active_set[g.a].append(g)

            for c in g.c:
                # get active pins on each control pin
                for anode in filter(None, active_set[c]):
                    sk.add_edge(anode, g, wt=1)
                    # graph[anode].add(g)

                if len(control_set[c]) > 0 and control_set[c][-1] is None:
                    control_set[c] = [g]
                else:
                    control_set[c].append(g)
                active_set[c].append(None)

        # sk = SkeletonGraph(graph, gates=ckt.gates, wires=w)
        sk.gates = ckt.gates
        sk.wires = w
        # handle edge weights

        # nx.set_edge_attributes(sk, 1, name='wt')

        # if augmented:
        #     # for each wire
        #     for u in range(w):
        #         # for each gate on this wire
        #         print(f"wire {u}: {sk.pred[u]}")
        #         for v in sk.pred[u]:
        #             # add an wire edge. positive for active, negative for
        #             # control. Remove regular edge.
        #             # del sk[v][u]['wt']
                    # sk[v][u]['wt'] = 1 if u == v.a else 1j

        # pprint(list(sk.edges(data=True)))
        return sk
    
    # def wire_canonical(self):
    #     topo = list(nx.topological_generations(self))

    #     n2l = {}

    #     for layer, nodes in enumerate(topo):
    #         n2l.update({n: layer for n in nodes})
        
    #     print(n2l)

    #     for layer, nodes in enumerate(topo):
    #         wire_nodes = list(filter(lambda t: type(t) == int, nodes))
    #         # concat all d['w']
    #         ie = defaultdict(list)
    #         for wn in wire_nodes:
    #             ie[wn] = tuple((n2l[u], u.type.value, self.out_degree(u, weight='weight'), e['weight']) for (u, v, e) in self.in_edges(wn, data=True))

    #         print(f"L{layer}: {ie}")
    #         assert(len(set(ie.values())) == len(ie))

    #         # print(f"{layer}: {list(filter(lambda t: type(t) == int, nodes))}")

    # TODO: need to differentiate between C-A collisions and A-C collisions.

    def label_in_deg(self, n, label):
        return sum(
            self[u][n].get(label, 0) for u in self.pred[n]
        )
    
    def sorted_gates(self):
        topo = list(nx.topological_generations(self))

        out = []
        for layer in topo:
            out.extend(sorted(layer))

        return out

    def canonical(self):
        topo = list(nx.topological_generations(self))

        n2l = {}
        l2n = defaultdict(list)
        types = {}

        for layer, nodes in enumerate(topo):
            n2l.update({n: layer for n in nodes})
            l2n[layer].extend(list(filter(lambda t: type(t) != int, nodes)))
            types.update({n: 'W' if type(n) == int else n.type.value for n in nodes})

        inD = self.in_degree(weight='wt')
        outD = self.out_degree(weight='wt')

        nx.set_node_attributes(self, None, name='hash')

        for layer, level in enumerate(topo):
            for node in level:
                mk = None
                if type(node) == int:
                    # get counts of node types per layer
                    mk = tuple(sum(g.a == node for g in lns) for l, lns in l2n.items())

                pred = tuple(sorted([self.nodes[u]['hash'] for u in self.pred[node]]))
                
                me = tuple((layer, types[node], inD[node], outD[node]))
                
                # print(f"{node} @ L{layer}:\n  {pred} {mk} {me}")
                h = hash((me, mk, pred))
                self.nodes[node]['hash'] = h

        # for n in self.nodes:
        #     # print(f"{n} {self.nodes[n]['hash']}")
        #     pass

        # for w in sorted(
        #     filter(lambda t: type(t) == int, self.nodes),
        #     key=lambda n: self.nodes[n]['hash']
        # ):
        # assert all_unique([self.nodes[n]['hash'] for n in filter(lambda t: type(t) == int and len(self.pred[t]) > 0, self.nodes)])

    def _colors(self):
        def _color_sel(n):
            if type(n) == int:
                # wires (augmented graph)
                return 'lightblue'
            else:
                if n.inverted:
                    return 'k'
                else:
                    return 'w'
        return list(map(_color_sel, self.nodes))

    def _label(self, g):
        if type(g) == int:
            return g
        else:
            return str(g).replace("[", "\n")[:-1]

    def r_neighborhood(self, v, r):
        ngh = nx.single_source_shortest_path_length(self, v, cutoff=r).keys()
        return list(ngh)

    def deg_profile(self):
        D = defaultdict(tuple)
        last = None
        for o in it.count():
            dp = self._deg_profile(o)
            dpv = sorted(dp.values())
            print(dpv)
            if dpv == last:
                print("no change, exit")
                break

            last = dpv

            for n, r in dp.items():
                if D[n]:
                    D[n] = tuple((*D[n], r))
                else:
                    D[n] = r

            if len(set(D.values())) == len(D):
                print("unique, exit")
                break

        print(f"done @ order {o}")
        print(D)
        return tuple(sorted(D.values()))

    def _deg_profile(self, order):
        """
        degree profile for isomorphism testing.
        order = 1: all nodes and their degrees
        order = 2: degrees of all neighbors

        for order r, this gives the degree profile of the r-neighborhood of each
        node.
        """

        inD = self.in_degree(weight="weight")
        outD = self.out_degree(weight="weight") 
 
        all_degrees = {n: (EdgeWeight(inD[n]), EdgeWeight(outD[n])) for n in self.nodes}

        # print(all_degrees)

        def node_type(n):
            return NodeType.Wire if type(n) == int else NodeType.Gate,

        # sort by manhattan distance first, then real, then imag (tuple)
        def norm(c):
            re = c.real
            im = c.imag
            return (abs(re) + abs(im), re, im)

        L = {}
        for v in self.nodes:
            ws = self.r_neighborhood(v, order)

            t = tuple(
                sorted(
                    [all_degrees[w] for w in ws],
                    key=lambda x: (node_type(x), norm(x[0]), norm(x[1])),
                )
            )
            L[v] = t

        return L

    def draw(self):
        nx.set_node_attributes(self, -10, name='layer')

        # sub = self.subgraph(self.graph['gates'])
        gen = list(nx.topological_generations(self))

        for layer, nodes in enumerate(gen):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                self.nodes[node]["layer"] = layer

        pos = nx.multipartite_layout(self, subset_key="layer")

        # xp = 0
        # min_x = min(map(lambda x: x[0], pos.values()))
        # max_x = max(map(lambda x: x[0], pos.values()))
        # min_y = min(map(lambda x: x[1], pos.values()))

        # space = (max_x - min_x) / self.graph['wires']

        # for n, (x, y) in pos.items():
        #     if type(n) == int:
        #         pos[n] = (xp, min_y)
        #         xp += space

        nx.draw(
            self,
            pos=pos,
            with_labels=False,
            connectionstyle="arc3,rad=0.4",
            node_color=self._colors(),
            edgecolors='k', # node outline
            edge_color=[random.choice(plt.cm.Dark2.colors) for _ in self.nodes],
            node_size=1200,
        )

        labels = {node: self._label(node) for node in self}
        nx.draw_networkx_labels(
            self,
            pos=pos,
            labels=labels,
            font_size=8,
            font_color='gray',
            font_weight="bold",
        )


class Permutation(tuple):
    def __init__(self, t):
        self.t = t

    def __str__(self):
        # TODO: cycle notation
        return (
            "["
            + " ".join(
                map(
                    lambda x: "." if x[0] == x[1] else str(x[1]),
                    zip(range(len(self.t)), self.t),
                )
            )
            + "]"
        )


def ckt_to_skel(n, ckt):
    return SkeletonGraph.from_circuit(Circuit(ckt, wires=n)).canonical()

def append_gates(base, n, skel):
    sk = list(it.chain.from_iterable(skel))
    r = []
    for b in base:
        ckt = Circuit(list((*sk, b)), wires=n)
        r.append(SkeletonGraph.from_circuit(ckt).canonical())
    return r

def perm_work(skel):
    # ignore the original circuit and recompute a canonical circuit from the
    # graph
    canon_ckt =  Circuit(list(it.chain.from_iterable(skel)))
    return skel, canon_ckt.perm()

def bitswaps(n):
    # skip first one, because that's the identity swap
    return list(it.permutations(range(n), n))[1:]

# TODO: is this right? this might only bitswap the input (or output, respectively)
# and we need to consider bitswap of both. e.g.
# ((N[0], N[1], N[2]), (T[2; 0 1],))
# and
# ((N[0], N[1], N[2]), (T[1; 0 2],))
# should be the same under bitswapping. I think we need to swap bits AND output
# positions.
def bitswap_iter(p, swaps):
    pa = np.array(p, dtype=int)
    
    for pi in swaps:
        out = np.zeros_like(pa, dtype=int)
        for i, b in enumerate(pi):
            out |= ((pa >> b) & 1) << i
        yield tuple(out.astype(list))

class SkeletonCache:
    PROCS = os.cpu_count() // 2

    def __init__(self, n, m):
        self.n = n
        self.m = m

        self.b = list(base_perms(n))
        # TODO: up to reordering
        self.n_perms = math.factorial(2**n) // 2
        self.num_ckts = len(self.b) ** m

    def build_from(self, prev):
        print(f"Building Skeleton Cache with {self.n=}, {self.m=} from scratch")
        assert prev.n == self.n and prev.m == self.m - 1

        self.num_ckts = len(prev.skel) * len(self.b)
        print(f"  {self.num_ckts} circuits")

        self.skel = Counter()
        self.perms = defaultdict(list)
        self.uniq_perms = defaultdict(list)

        with Progress(
            "[progress.description]{task.description}",
            MofNCompleteColumn(),
            BarColumn(),
            "[progress.percentrage]{task.percentage:>3.0f}%",
            TimeElapsedColumn(),
            "ETA",
            TimeRemainingColumn(),
            TransferSpeedColumn(),
            refresh_per_second=4,
            speed_estimate_period=120
        ) as prog:
            
            task_skel = prog.add_task("Skeletons", total=self.num_ckts)

            with mp.Pool(self.PROCS) as pool:
                par_it = pool.imap_unordered(
                    partial(append_gates, self.b, self.n), prev.skel.keys()
                )

                for c in par_it:
                    prog.advance(task_skel, advance=len(c))
                    self.skel.update(c)

            task_perm = prog.add_task("Perms", total=len(self.skel))
            with mp.Pool(self.PROCS) as pool:
                par_it = pool.imap_unordered(
                    perm_work, self.skel.keys(), chunksize=16
                )

                # TODO: inverse function can be ignored
                for skel, perm in par_it:
                    prog.advance(task_perm)
                    self.perms[perm].append(skel)

            swaps = bitswaps(self.n)
            
            if self.n > 5:
                print(f"Skipping: {self.n=} with {len(swaps)} bit swaps")
                self.uniq_perms = {}
            else:
                # single proc for now
                task_uniq_perm = prog.add_task("Canon. Perm", total=len(self.perms))
                for perm, skels in self.perms.items():
                    prog.advance(task_uniq_perm)
                    pp = perm
                    for p in bitswap_iter(perm, swaps):
                        # print(f" b/s: {p}")
                        if p in self.uniq_perms:
                            pp = p
                            # print(f"duplicate: {perm} => {p} already seen")
                            break

                    # TODO: recanonicalize these skeletons?
                    self.uniq_perms[pp].extend(skels)
            
            prog.refresh()

    def build(self):
        print(f"Building Skeleton Cache with {self.n=}, {self.m=} from scratch")

        print(
            f"2^{round(math.log2(self.n_perms))} perms on {self.n} wires. b={len(self.b)}; {self.num_ckts} circuits."
        )
        skel = Counter()
        perms = defaultdict(list)
        uniq_perms = defaultdict(list)

        with Progress(
            "[progress.description]{task.description}",
            MofNCompleteColumn(),
            BarColumn(),
            "[progress.percentrage]{task.percentage:>3.0f}%",
            TimeElapsedColumn(),
            "ETA",
            TimeRemainingColumn(),
            TransferSpeedColumn(),
            refresh_per_second=4,
            speed_estimate_period=120
        ) as prog:
            task_skel = prog.add_task("Skeleton Graphs", total=self.num_ckts)

            with mp.Pool(self.PROCS) as pool:
                par_it = pool.imap_unordered(
                    partial(ckt_to_skel, self.n), all_circuits(self.n, self.m), chunksize=16
                )

                for canon in par_it:
                    prog.advance(task_skel)

                    # For distribution, do we care about the number?
                    # Sample uniformly from _all skeletons_ or _all circuits_
                    skel[canon] += 1

                self.skel = skel

            task_perm = prog.add_task("Permutations", total=len(skel))
            with mp.Pool(self.PROCS) as pool:
                par_it = pool.imap_unordered(
                    perm_work, skel.keys(), chunksize=16
                )

                # TODO: inverse function can be ignored
                for skel, perm in par_it:
                    prog.advance(task_perm)

                    perms[perm].append(skel)

                self.perms = perms

            swaps = bitswaps(self.n)
            
            if self.n > 5:
                print(f"Skipping: {self.n=} with {len(swaps)} bit swaps")
                self.uniq_perms = {}
            else:
                # single proc for now
                task_uniq_perm = prog.add_task("Canon. Perm", total=len(perms))
                for perm, skels in perms.items():
                    prog.advance(task_uniq_perm)
                    pp = perm
                    for p in bitswap_iter(perm, swaps):
                        # print(f" b/s: {p}")
                        if p in uniq_perms:
                            pp = p
                            # print(f"duplicate: {perm} => {p} already seen")
                            break

                    # TODO: recanonicalize these skeletons?
                    uniq_perms[pp].extend(skels)

                self.uniq_perms = uniq_perms

            prog.refresh()

    def stats(self):
        # print(self.skel)
        # print(self.perms)
        num_graphs = len(self.skel)
        num_perms = len(self.perms)

        # iden_perm = Permutation(range(2**self.n))
        # iden_graphs = self.G[iden_perm]
        # num_id_graphs = len(iden_graphs)
        # num_id_ckts = sum(iden_graphs.values())

        # uniq_graph = sum(1 if len(v) == 1 else 0 for v in self.G.values())
        # uniq_ckt = sum(1 if (sum(v.values()) == 1) else 0 for v in self.G.values())

        perm_pct = num_perms / self.num_ckts * 100

        # stored = sum(sum(g.values()) for g in self.G.values())

        print(
            f"{self.num_ckts} circuits => {num_graphs} graphs => {num_perms} perms ({perm_pct:.1f}%)"
        )

        dups = num_perms - len(self.uniq_perms)
        if num_perms > 0:
            dup_pct = 100 * dups / num_perms
            print(f"  {len(self.uniq_perms)} uniq perms up to wire ordering ({dups=}, {dup_pct:.1f}%)")
        
        # # print(f"  Stored {stored}")
        # print(f"  {num_id_ckts} id ckts => {num_id_graphs} id graphs")
        # print(f"  {uniq_graph} perms w/ a uniq graph")
        # print(f"  {uniq_ckt} perms w/ a uniq circuit")

        n_id = 0
        for n_ in range(1, self.n):
            pow = 2 ** n_
            p = tuple(range(pow))
            # print(f"{p}: {len(self.uniq_perms[p])} skel")
            n_id += len(self.uniq_perms.get(p, []))

        print(f"  {n_id} identity graphs")

        # ok, well up to reordering this is not super accurate. just a good
        # heuristic
        n_uniq = sum(1 if len(v) == 1 else 0 for v in self.uniq_perms.values())
        print(f"  {n_uniq} perms with a unique skeleton")

        # for p, sk in self.uniq_perms.items():
        #     if len(sk) == 1:
        #         print(f"Perm {p}")
        #         skel = sk[0]
        #         # print(skel)
        #         # print(Circuit(list(it.chain.from_iterable(skel))))

    def print(self):
        sG = sorted(self.G.items(), key=lambda x: sum(x[1].values()), reverse=True)
        for p, gs in sG:
            title = f"Perm {p}: {len(gs)} graphs, {sum(gs.values())} circuits"
            print("\n" + "=" * len(title))
            print(title)

            for i, (g, ckt) in enumerate(gs.items()):
                print(f"#{1+i} Graph {g}: {ckt} circuits.")

                print(Circuit(list(it.chain.from_iterable(g)), wires=self.n))

                for k in ckt:
                    print(k)
                    break
