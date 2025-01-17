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
from rich import print
import multiprocessing
import os

class GateT(Enum):
    NOT = "N"
    CNOT = "C"
    CCNOT = "T"


class Gate:
    """A general CxNOT gate."""

    c = []

    def __init__(self, a, c1=None, c2=None):
        assert type(a) is int
        if c1 is None:
            # NOT gate
            self.a = a
            self.c = ()
            self.type = GateT.NOT
        elif c2 is None:
            # CNOT gate
            assert type(c1) is int
            self.a = a
            self.c = (c1,)
            self.type = GateT.CNOT
        else:
            # CCNOT gate
            assert type(c1) is int and type(c2) is int
            self.a = a
            self.c = (c1, c2)
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
        name = self.type.value if not show_id else f"{self.type.value}#{id(self)}"
        if self.type == GateT.NOT:
            return f"{self.type.value}[{self.a}]"
        return f"{name}[{self.a}; {' '.join(map(str, self.c))}]"

    def compact(self):
        short_type = {GateT.NOT: "N", GateT.CNOT: "C", GateT.CCNOT: "T"}[self.type]
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


class CanonicalGate(Gate):
    # def __init__(self, a, c1=None, c2=None):
    #     super().__init__(a, c1, c2)

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
        self.gates.extend(copy.deepcopy(g) for g in gates)
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
        return "{" + " ".join(g.compact() for g in self.gates) + "}"


def all_circuits(n, m):
    b = base_perms(n)
    yield from it.product(b, repeat=m)

class PermutationDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._next_value = 0

    def __missing__(self, k):
        while self._next_value in self.values():
            self._next_value += 1
        self[k] = self._next_value
        return self._next_value


class SkeletonGraph(nx.DiGraph):
    """A lightweight class to represent skeleton graphs."""

    # edges do not need labels
    edge_dict = {}

    def from_circuit(ckt: Circuit):
        w = ckt.wires
        active_set = [[] for _ in range(w)]
        control_set = [[] for _ in range(w)]

        graph = {g: set() for g in ckt.gates}

        for g in ckt.gates:
            # get control pins on this active pin
            back = control_set[g.a]
            for h in filter(None, back):
                graph[h].add(g)

            control_set[g.a].append(None)

            if len(active_set[g.a]) > 0 and active_set[g.a][-1] is None:
                active_set[g.a] = [g]
            else:
                active_set[g.a].append(g)

            for c in g.c:
                # get active pins on each control pin
                for anode in filter(None, active_set[c]):
                    graph[anode].add(g)

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

        new_topo = []

        for level in topo:
            ls = sorted(level)
            t = []
            for g in ls:
                t.append(CanonicalGate(*[wire_perm[p] for p in g.wires()]))

            new_topo.append(tuple(sorted(t)))

        # pprint(wire_perm)
        return tuple(new_topo)

    def draw(self):
        gen = list(nx.topological_generations(self))

        for layer, nodes in enumerate(gen):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                self.nodes[node]["layer"] = layer

        pos = nx.multipartite_layout(self, subset_key="layer")

        nx.draw(
            self,
            pos=pos,
            with_labels=False,
            connectionstyle="arc3,rad=0.4",
            node_color="black",
            edge_color=[random.choice(plt.cm.Dark2.colors) for _ in self.nodes],
            node_size=1200,
        )

        labels = {node: str(node).replace("[", "\n")[:-1] for node in self}
        nx.draw_networkx_labels(
            self,
            pos=pos,
            labels=labels,
            font_size=8,
            font_color="white",
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


    def build(self):
        print(f"Building Skeleton Cache with {self.n=}, {self.m=}")
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

            with multiprocessing.Pool(self.PROCS) as pool:
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
            with multiprocessing.Pool(self.PROCS) as pool:
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
                print(f"Skipping: {self.n=} and {swaps=}")
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
            n_id += len(self.uniq_perms[p])

        print(f"  {n_id} identity graphs")

        # ok, well up to reordering this is not super accurate. just a good
        # heuristic
        n_uniq = sum(1 if len(v) == 1 else 0 for v in self.uniq_perms.values())
        print(f"  {n_uniq} perms with a unique skeleton")

        # for p, sk in self.uniq_perms.items():
        #     if len(sk) == 1:
        #         print(f"Perm {p}")
        #         skel = sk[0]
        #         print(skel)
        #         print(Circuit(list(it.chain.from_iterable(skel))))

    def print(self):
        sG = sorted(self.G.items(), key=lambda x: sum(x[1].values()), reverse=True)
        for p, gs in sG:
            title = f"Perm {p}: {len(gs)} graphs, {sum(gs.values())} circuits"
            print("\n" + "=" * len(title))
            print(title)

            for i, (g, ckt) in enumerate(gs.items()):
                print(f"#{1+i} Graph {g}: {ckt} circuits.")

                print(Circuit(list(it.chain.from_iterable(g)), wires=self.n))

                # for k in ckt:
                #     print(k)
                # break
