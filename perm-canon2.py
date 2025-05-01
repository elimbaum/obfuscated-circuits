#!/usr/bin/env python3
#
# A new method for canonicalizing permutations according to lexicographical
# ordering. However, this only requires evaluating the permutation (circuit)
# at O(n), rather than O(n!), locations, saving (n-1)!x the work.
#
# To keep the algorithm simple, we define a new ordering of the _indices_ of the
# permutation which makes canonicalization easier. Consider the indices of the
# permutation according to how "free" their bits are to move; that is, order
# the permutation indices i according to (2n choose wt(i)), where wt(i) returns
# the hamming weight of i. This change drastically simplifies the algorithm for
# reordering permutation indices.
#
# Concretely, we consider permutation indices
# 0, (2^n - 1)     [no bit shuffling possible]
# (2^i, ...)       [single 1 bit]
# (2^n - 2^i, ...) [single 0 bit]
# etc...           [two 1 bits]
#
# In this way, it is extremely unlikely that we will need to consider
# permutation indices with a hamming weight of 2 or more in order to fully
# canonicalize a permutation.

import os
import sys
import time
import numpy as np
from random import shuffle, seed, randint
from rich import print
import pprint
from functools import reduce
import operator
import copy
from collections import defaultdict, Counter

s_ = randint(0, 1e6)
s_ = 169677
print("SEED:", s_)
seed(s_)

def strings_of_weight(w):
    # due to Bill Gosper?
    def next(x):
        c = x & -x
        r = x + c
        return (((r ^ x) >> 2) // c) | r

    # w 1's in the LSB
    alpha = int(2**w) - 1

    if w == 0:
        yield alpha
        return

    while alpha < N:
        yield alpha
        alpha = next(alpha)


def index_sets():
    for s in range(n // 2 + 1):
        # "light" string (more 0s than 1s)
        yield strings_of_weight(s)

        if 2 * s == n:
            continue

        # "heavy" string (bitwise complement, more 1s than 0s)
        yield strings_of_weight(n - s)


def apply_shuffle(p, sh):
    k = np.zeros_like(p, dtype=int)
    i1 = np.arange(len(p))
    i2 = np.zeros_like(i1)

    for (
        src,
        dest,
    ) in enumerate(sh):
        k |= ((p >> src) & 1) << dest
        i2 |= ((i1 >> src) & 1) << dest

    k = k[i2.argsort()]
    return k


class CandSet:
    def __init__(self, n):
        self.map = {i: set(range(n)) for i in range(n)}
        self.n = n

    def preimages(self, w):
        # what are the bits of w?
        wbits = set(i for i in range(self.n) if (w & (1 << i)))
        # print(f"{wbits=}")

        # for which equal-weight s does there still exist a valid candidate
        # mapping from s to w?
        p = set()
        for s in strings_of_weight(w.bit_count()):
            # source bits map to target wbits
            sbits = [i for i in range(self.n) if (s & (1 << i))]

            # may not be totally correct, just a first pass
            # however might be correct if all distinct mappings are disjoint...
            # tbd
            bit_possib = reduce(operator.or_, (self.map[_b] for _b in sbits), set())

            if wbits <= bit_possib:
                p.add(s)

        return p

    def __str__(self):
        return pprint.pformat(self.map)

    def normalize(self):
        changes = False
        for i in range(n):
            c = [k for k, v in self.map.items() if i in v]
            # if my only assignment, no one else can have it
            if len(c) == 1:
                changes = True
                self.map[c[0]] = {i}

            my_matches = self.map[i]
            dups = sum(my_matches == v for v in self.map.values())

            if len(my_matches) == dups:
                changes = True
                for k, v in self.map.items():
                    if v == my_matches:
                        continue
                    v -= self.map[i]

        return changes

    def complete(self):
        return all(len(v) == 1 for v in self.map.values())

    def legal(self):
        return all(len(v) > 0 for v in self.map.values())

    def mapto(self, i):
        return set(k for k, v in self.map.items() if i in v)

    def fixmap(self, from_, to):
        for k, v in self.map.items():
            if k == from_:
                v &= {to}
            else:
                v.discard(to)

    def isect(self, other):
        for k, v in self.map.items():
            v &= other.map[k]

    def union(self, other):
        for k, v in self.map.items():
            v |= other.map[k]

    def enforce(self, x, y):
        # enforce x -> y
        x_zeros = set(j for j in range(n) if (x & (1 << j)) == 0)
        x_ones = set(j for j in range(n) if (x & (1 << j)) != 0)
        y_zeros = set(j for j in range(n) if (y & (1 << j)) == 0)
        y_ones = set(j for j in range(n) if (y & (1 << j)) != 0)

        for k, v in self.map.items():
            if k in x_zeros:
                v &= y_zeros

            if k in x_ones:
                v &= y_ones

    def output(self):
        if self.complete():
            return [self.map[i].pop() for i in range(n)]
    
        # incomplete. any assignment is valid; choose one at random
        # TODO: is this ok? only appears to happen with n=3
        print("Incomplete output.")
        L = []
        for i in range(n):
            p = min(self.map[i])
            L.append(p)
            self.fixmap(i, p)
            print(self.map)

        return L

    # does the mapping allow x -> y?
    def is_consistent(self, x, y):
        # print(f"Check consistent: {x} -> {y}")
        x_zeros = set(j for j in range(n) if (x & (1 << j)) == 0)
        x_ones = set(j for j in range(n) if (x & (1 << j)) != 0)
        y_zeros = set(j for j in range(n) if (y & (1 << j)) == 0)
        y_ones = set(j for j in range(n) if (y & (1 << j)) != 0)

        # print(f"{x_zeros=}, {y_zeros=}")
        # print(f"{x_ones=}, {y_ones=}")

        # print("is_con:", self.map)

        for k, v in self.map.items():
            # mapping from a zero... to a zero?
            if k in x_zeros:
                if len(y_zeros & v) == 0:
                    return False

            # mapping from a one to a one?
            if k in x_ones:
                if len(y_ones & v) == 0:
                    return False

        return True
    
    def __lt__(self, other):
        for i in range(self.n):
            s = self.map[i]
            o = other.map[i]

            if len(s) < len(o):
                return True
            elif len(s) > len(o):
                return False
            # same length.
            elif self.map[i] < other.map[i]:
                return True
            elif self.map[i] > other.map[i]:
                return False
            
        return False


# maps INPUT BITS to POSSIBLE OUTPUT BITS
# we terminate when each row (input bit) has a single non-null output


def min_consistent(x, cand):
    # zeros in the input
    zlist = set(j for j in range(n) if (x & (1 << j)) == 0)
    # ones
    olist = set(range(n)) - zlist

    c2 = copy.deepcopy(cand)

    a = {}
    for i in range(n - 1, 0 - 1, -1):
        # which inputs bits could map to this one?
        i_from = c2.mapto(i)

        if len(i_from) == 0:
            # print(f"{i=} -> {i_from=}")
            # print(c2)
            return None, None

        # are any of them zeros?
        cand_z = i_from & zlist
        cand_o = i_from & olist
        # print(f"{i=} {i_from=}")
        # j = None

        if len(cand_z):
            # print("Got zero")
            # most significant zero
            j = max(cand_z)
            a[j] = i
        else:
            # least significant one
            j = min(cand_o)

        a[j] = i
        # print(f"sending a[{j}] to {i}")
        c2.fixmap(j, i)
        # print(c2)

    # print(c2)

    out = 0
    for k, v in a.items():
        # print("bit", k, v)
        out |= ((x >> k) & 1) << v
        # print(bin(out))

    # print(f"Try to map {x} -> {out} b{out:0{n}b}")

    zout = set(j for j in range(n) if (out & (1 << j)) == 0)

    c3 = CandSet(n)
    for k, v in c3.map.items():
        if k in zlist:
            v &= zout
        else:
            v -= zout

    # print("c3", c3)

    return out, c3

end_w = None

def canonicalize(perm, do_print=True):

    prnt = print
    if not do_print:
        prnt = lambda x, *args: 0

    n = (len(perm) - 1).bit_length()
    cand = CandSet(n)

    prnt("Input:", perm)

    steps = 0

    for s in index_sets():
        # list of perm indices with a given weight, according to our defined order
        sl = list(s)

        prnt(f"===========================")
        prnt(f"==== Index set: {sl}")
        prnt(f"===========================")

        # for each perm index in this weight class
        for w in sl:
            steps += 1
            if n > 3 and (w.bit_count() == 2 or w.bit_count() == (n - 2)):
                prnt("Skip multibit candidates. Aborting.")
                # os.abort()

            # what are the valid preimages of w?
            prnt("Candidate mappings:")
            prnt(cand)

            p = cand.preimages(w)

            prnt(f"\n--- π({w}) = {perm[w]} b{perm[w]:0{n}b}")
            prnt(f"    preimage({w}) = {p}")

            if len(p) > 1:
                for x in p:
                    if x == w:
                        continue
                    prnt(f"     π({x}) = {perm[x]} b{perm[x]:0{n}b}")

            passed = []
            best = None
            best_x = None
            identity = False
            for x in p:
                y = perm[x]
                cand2 = copy.deepcopy(cand)
                cand2.enforce(x, w)
                cand2.normalize()
                val, m = min_consistent(y, cand2)

                if val == None:
                    continue

                # print(f"p[{x}] = {y}. Consistency map for min={val}:")
                # print(m)
                m.isect(cand)
                # print("Isect.")
                # print(m)

                # check if valid
                is_con = m.is_consistent(x, w)
                prnt(f"  Consistent π({x}) -> π'({w})? {is_con}")

                if not is_con:
                    continue

                prnt(f"π({x=})=(y={int(y)}); {w=} val={int(val)}")

                if best is None or val < best:
                    if best is not None:
                        prnt(f"    {y} -> {val} < {best}")
                    else:
                        prnt(f"    {y} -> {val}")
                    best = val
                    best_x = x
                    passed = [m]
                    if w == val:
                        identity = True
                elif val == best:
                    if w == val:
                        prnt(f"    {y} -> {val} == {best}, replacing identity")
                        if identity:
                            passed.append(m)
                        else:
                            passed = [m]

                        identity = True
                        best_x = x
                    elif not identity:
                        prnt(f"    {y} -> {val} == {best}, adding")
                        passed.append(m)
                else:
                    prnt(f"    {y} -> {val} > {best}, skip")
                    continue

            if len(passed) == 0:
                # no valid candidate passed checks; skip
                continue

            # if multiple matches, no unique lex assignment. move on to next
            # word... or do union?
            if len(passed) > 1:
                prnt("Warning: ambiguous! Choosing from")
                for p in passed:
                    prnt(f"p={p}")

                # continue
                m = passed[0]
                for p in passed[1:]:
                    m.union(p)

                prnt("union=", m)
                ## These lines fix the hard-identity case but break the random
                ## case
                # cand.isect(m)
                # cand.normalize()
            else:
                prnt(f"  Updating candidate map: π({best_x}) -> π'({w})")
                cand = passed[0]

            prnt("---")

            if cand.complete():
                break

        if cand.complete():
            break

    prnt("Final map:")
    prnt(cand)

    # if we made it outside the loop without canonicalizing, abort
    if not cand.complete():
        prnt("Incomplete!")
        
        # os.abort()

    global end_w
    end_w = steps

    return cand.output()


n = 5
R = np.arange(n)
N = 2**n
bincomp = 1 << R

perm = np.arange(N)
shuffle(perm)

s = list(strings_of_weight(2))

# perm[s[-1]], perm[s[-2]] = perm[s[-2]], perm[s[-1]]
# perm[s[-3]], perm[s[-4]] = perm[s[-4]], perm[s[-3]]

print(perm)

print("Perm:", perm)

q = canonicalize(perm)
pp = apply_shuffle(perm, q)
print(perm, q, pp)
assert set(q) == set(range(n))
assert set(pp) == set(range(N))

print(end_w)

test = True

times = []

if test:
    good = 0
    max_w = Counter()
    for _ in range(10000):
        print(_)
        shuffle(perm)
        start = time.time()
        q = canonicalize(perm, do_print=False)
        end = time.time()
        times.append(end-start)

        # print(q)
        assert set(q) == set(range(n))
        assert set(pp) == set(range(N))

        print(f"{q=}")

        pp = apply_shuffle(perm, q)
        max_w[end_w] += 1
        for _ in range(4):
            shuffle(q)
            p2 = apply_shuffle(perm, q)

            if np.array_equal(perm, p2):
                continue

            q2 = canonicalize(p2, False)
            pp2 = apply_shuffle(p2, q2)

            good += 1

            # print("RESULT")
            # print("expect", pp)
            # print("got   ", pp2)
            assert np.array_equal(pp, pp2)
            # print("## CORRECT ##")

        print(f"{max_w=}")

    print(f"{good=}")

    print(np.mean(times) * 1e6)
