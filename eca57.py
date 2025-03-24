#!/usr/bin/env python3
from functools import reduce
from rich import print
from rich.progress import track
import readline
import itertools as it
import random
from collections import Counter
import math

class ECA():
    def __init__(self, a, c1, c2):
        self.c = (c1, c2)
        self.a = a

    def top(self):
        return min(self.a, min(self.c))

    def bottom(self):
        return max(self.a, max(self.c))
    
    def eval(self, av, c1v, c2v):
        return av ^ (c1v or (not c2v))

class Circuit():
    def __init__(self, gates):
        self.gates = gates
        self.wires = 1 + max(g.bottom() for g in gates)

    def from_string(s):
        gate_str = s.split(";")
        gates = []
        for g in gate_str:
            g = g.strip()
            if not g:
                continue
            a, c1, c2 = tuple(map(int, g.split()))
            gates.append(ECA(a, c1, c2))

        return Circuit(gates)
    
    def perm(self, n=None):
        p = []
        if n is None:
            n = self.wires

        for i in range(2**n):
            b = format(i, "0{}b".format(n))
            out_list = self.eval(list(map(int, b)))

            out_n = reduce(lambda a, b: (a << 1) | b, out_list, 0)
            p.append(out_n)

        return tuple(p)
    
    def eval(self, input):
        # assert len(input) <= self.wires

        state = input[:]
        for g in self.gates:
            state[g.a] = g.eval(state[g.a], *(state[c] for c in g.c))

        return state

    def __str__(self):
        CKT_WIRE_CHAR = "─"
        CKT_CONTROL_CHAR = "│"
        CKT_DOT_CHAR = "●"
        CKT_OPEN_DOT = "○"
        s = ""
        for w in range(self.wires):
            s += f"{w:>3} {CKT_WIRE_CHAR}"
            for g in self.gates:
                if w == g.a:
                    s += "( )"
                elif w == g.c[0]:
                    s += CKT_WIRE_CHAR + CKT_DOT_CHAR + CKT_WIRE_CHAR
                elif w == g.c[1]:
                    s += CKT_WIRE_CHAR + CKT_OPEN_DOT + CKT_WIRE_CHAR
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
        return ";".join(f"{g.a} {g.c[0]} {g.c[1]}" for g in self.gates)
    
    def __len__(self):
        return len(self.gates)

def no_pairwise_adj(B, m):
    for g in it.product(B, repeat=m):
        if any(x == y for x, y in it.pairwise(g)):
            # some g_i == g_{i+1}, so skip this circuit
            continue

        yield g

class CircuitRecord():
    ckt = None
    count = 0

    def __init__(self, name):
        self.name = name

    def add(self, ckt):
        if self.ckt is None:
            self.ckt = ckt
            self.count = 1
        elif len(ckt) < len(self.ckt):
            self.ckt = ckt
            self.count = 1
        elif len(ckt) == len(self.ckt):
            self.count += 1
    
    def __str__(self):
        return f"{self.name} ({self.count} circuits)"


special_perms = {
    (0, 1, 2, 3, 4, 5, 6, 7): CircuitRecord("Identity"),
    (0, 1, 2, 3, 4, 5, 7, 6): CircuitRecord("AND"),
    (1, 0, 3, 2, 5, 4, 6, 7): CircuitRecord("NAND"),
    (0, 1, 3, 2, 5, 4, 7, 6): CircuitRecord("OR"),
    (1, 0, 2, 3, 4, 5, 6, 7): CircuitRecord("NOR"),
    # (a & !b)
    (0, 1, 2, 3, 5, 4, 6, 7): CircuitRecord("Implies"),
    # (a & !b)
    (1, 0, 3, 2, 4, 5, 7, 6): CircuitRecord("NImplies"),

    (0, 1, 3, 2, 5, 4, 6, 7): CircuitRecord("XOR"),
    (1, 0, 2, 3, 4, 5, 7, 6): CircuitRecord("XNOR"),
    (1, 0, 3, 2, 5, 4, 7, 6): CircuitRecord("NOT"),
    # (a)
    (0, 1, 2, 3, 5, 4, 7, 6): CircuitRecord("CNOT"),
    # (a)
    (1, 0, 3, 2, 4, 5, 6, 7): CircuitRecord("NCNOT"),

    (7, 6, 5, 4, 3, 2, 1, 0): CircuitRecord("Flip"),
}

WIRES = 3

def base_perms(n):
    return list(it.permutations(range(n), r=3))

def find_ckt(m, B=None):
    if not B:
        B = base_perms(WIRES)
    bb = len(B)

    perms = Counter()

    total = bb * (bb - 1) ** (m - 1)

    print(f"{m} gates, {total} circuits")

    for gl in track(no_pairwise_adj(B, m), total=total):
        ckt = Circuit([ECA(*g) for g in gl])
        p = ckt.perm(WIRES)
        # if p == tuple(range(8)):
        #     print(ckt)

        if p in special_perms:
            special_perms[p].add(ckt)

        perms[p] += 1

    print(perms.most_common(16)[::-1])
    
    pct = len(perms) / total * 100
    ppct = len(perms) / (math.factorial(2 ** WIRES) // 2) * 100
    print(f"{len(perms)} distinct perms ({pct:.1f}% of circuits, {ppct:.1f}% of perms)")

    return perms

q = None
last_ckt = None

def offset(a, c1, c2):
    if c1 > a and c1 > c2:
        return max(a, c2) - c1
    elif c1 < c2:
        return min(a, c2) - c1
    else:
        return a - c1

def cmd(q):
    global last_ckt
    if q.startswith("check"):
        m = int(q.split()[1])
        if m < 0:
            perms = set()

            for mm in range(1, -m + 1):
                p_ = find_ckt(mm).keys()
                if mm == -m:
                    diff = p_ - perms
                    print(diff)
                    print(f"{len(diff)} new perms")
                    
                perms.update(p_)
                print(f"{len(perms)=}")
                print()
        else:
            find_ckt(m)
    elif q.startswith("special"):
        for p, ckt in special_perms.items():
            print(f"{ckt} {p}")
            print(ckt.ckt)
            print(repr(ckt.ckt))
    elif q.startswith("wires"):
        global WIRES
        split = q.split()
        if len(split) > 1:
            WIRES = int(split[1])
        print(f"n={WIRES}")
    elif q.startswith("latex"):
        if not last_ckt:
            return

        print(last_ckt)

        W = last_ckt.wires
        wire_store = [[] for _ in range(W)]

        for g in last_ckt.gates:
            a = g.a
            c1, c2 = g.c

            for w in range(W):

                if w == a:
                    wire_store[w].append("\\targ")
                elif w == c1:
                    d = offset(a, c1, c2)
                    wire_store[w].append(f"\\ctrl{{{d}}}")
                elif w == c2:
                    d = offset(a, c2, c1)
                    wire_store[w].append(f"\\ctrlo{{{d}}}")
                else:
                    wire_store[w].append("\\qw")

        for w in range(W):
            s = " & ".join(wire_store[w])
            print(f"& {s} & \\qw \\\\")

    else:
        try:
            ckt = Circuit.from_string(q)
        except ValueError:
            print("parse error")
            return

        print()
        print(ckt)
        last_ckt = ckt
        p = ckt.perm()
        for i, j in enumerate(p):
            diffs = j ^ i
            os = f"{diffs:>03b}".replace("0", ".")
            print(f"{i:2}  {j:2} {j:>03b} {os}")

while not q or not q.startswith("q"):
    if q:
        cmd(q)
    
    try:
        q = input("> ").strip().lower()
    except EOFError:
        break
    except KeyboardInterrupt:
        q = None
        print()
        continue

