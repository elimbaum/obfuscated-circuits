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
import numpy as np
from random import shuffle, seed, randint
from rich import print
import pprint
from functools import reduce
import operator
import copy
from collections import defaultdict

s_ = randint(0, 1e6)
# s_ = 742732
# s_ = 755845
# s_ = 762105 # makes it to multibit
# s_ = 829866 # makes it to multibit
# s_ = 599345
print("SEED:", s_)
seed(s_)

n = 12
R = np.arange(n)
N = 2 ** n
bincomp = 1 << R

perm = np.arange(N)
shuffle(perm)

print("Perm:", perm)

def strings_of_weight(w):    
    # due to Bill Gosper?
    def next(x):
        c = x & -x
        r = x + c
        return (((r ^ x) >> 2) // c) | r

    # w 1's in the LSB
    alpha = int(2 ** w) - 1

    if w == 0:
        yield alpha
        return

    while alpha < N:
        # print(f"{alpha:7} {alpha:07b}")
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

    for src, dest, in enumerate(sh):
        k |= ((p >> src) & 1) << dest
        i2 |= ((i1 >> src) & 1) << dest

    k = k[i2.argsort()]
    return k

class CandidateMap():
    def __init__(self, n):
        self.a = np.ones(dtype=bool, shape=(n, n))
        self.n = n

    def __str__(self):
        s = ""
        spacing = int(np.ceil(np.log10(n)))
        for j, r in enumerate(self.a):
            s += f"{j:{spacing}} | "
            # q = r * np.arange(self.n)[::-1]
            s += ' '.join(f"{(n - i - 1 if qi else ''):{spacing}}" for i, qi in enumerate(reversed(r)))
            s += '\n'

        return s
    
    def normalize(self):
        # does this actually need 2 passes? more?
        for _ in range(2):
            row_sum = np.sum(self.a, axis=1)
            for r in np.where(row_sum == 1)[0]:
                d = np.where(self.a[r, :])[0]
                self.a[:, d] = False
                self.a[r, d] = True

            col_sum = np.sum(self.a, axis=0)
            for c in np.where(col_sum == 1)[0]:
                d = np.where(self.a[:, c])[0]
                self.a[d, :] = False
                self.a[d, c] = True


    def complete(self):
        row_sum = np.sum(self.a, axis=1)
        return all(row_sum == 1)
    
    def legal(self):
        col_sum = np.sum(self.a, axis=0)
        row_sum = np.sum(self.a, axis=1)

        return all(col_sum > 0) and all(row_sum > 0)

    
class CandSet():
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

    # does the mapping allow x -> y?
    def is_consistent(self, x, y):
        print(f"Check consistent: {x} -> {y}")
        x_zeros = set(j for j in range(n) if (x & (1 << j)) == 0)
        x_ones = set(j for j in range(n) if (x & (1 << j)) != 0)
        y_zeros = set(j for j in range(n) if (y & (1 << j)) == 0)
        y_ones = set(j for j in range(n) if (y & (1 << j)) != 0)

        print(f"{x_zeros=}, {y_zeros=}")
        print(f"{x_ones=}, {y_ones=}")

        for k, v in self.map.items():
            # mapping from a zero... to a zero?
            # print(f"Check {k}: {v}")
            if k in x_zeros:
                if len(y_zeros & v) == 0:
                    return False
                
            # mapping from a one to a one?
            if k in x_ones:
                if len(y_ones & v) == 0:
                    return False
            
        return True

# maps INPUT BITS to POSSIBLE OUTPUT BITS
# we terminate when each row (input bit) has a single non-null output

# cand = CandidateMap(n)

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

    print(f"Try to map {x} -> {out} b{out:0{n}b}")

    zout = set(j for j in range(n) if (out & (1 << j)) == 0)

    c3 = CandSet(n)
    for k, v in c3.map.items():
        if k in zlist:
            v &= zout
        else:
            v -= zout

    print("c3", c3)

    return out, c3

def canonicalize(p):
    n = (len(p) - 1).bit_length()
    cand = CandSet(n)

    print("Starting map")
    print(cand)

    # print(cand)

    for s in index_sets():
        # list of perm indices with a given weight, according to our defined order
        sl = list(s)

        print(f"==== Index set: {sl}")
        
        # if 15 in sl:
        #     print("skip single zero set")
        #     break

        if 3 in sl:
            print("skip multibit")
            os.abort()
            break

        for w in sl:
            # for each perm index in this weight class
            p = cand.preimages(w)

            print("============")
            print("============")
            print(f"π({w}) = {perm[w]} b{perm[w]:0{n}b}")
            print(f"pre({w}) = {p}")

            if len(p) > 1:
                print("---- others")
                for x in p:
                    if x == w:
                        continue
                    print(f"  π({x}) = {perm[x]} b{perm[x]:0{n}b}")
                print("----")

            
            passed = set()
            best = None
            for x in p:
                y = perm[x]
                val, m = min_consistent(y, cand)

                print(f"p[{x}] = {y}. Consistency map:")
                print(m)
                m.isect(cand)
                print("Isect.")
                print(m)

                # check if valid
                is_con = m.is_consistent(x, w)
                print(f"Consistent π({x}) -> π'({w})? {is_con}")

                if not is_con:
                    continue

                if best is None or val <= best:
                    best = val
                else:
                    print(f"{val} > {best}")
                    continue

                passed = {m}

            if len(passed) == 0:
                continue
            if len(passed) == 1:
                cand = passed.pop()
            else:
                print("Multiple:")
                for p in passed:
                    print(p)
                os.abort()

            print(cand)

            if cand.complete():
                break
            
        if cand.complete():
            break
    
    if not cand.complete():
        os.abort()

    return [cand.map[i].pop() for i in range(n)]


            # loop from MSB to LSB
            # removed = 0
    #         skipped = 0
    #         for i in range(n - 1, 0 - 1, -1):
    #             print("########")
    #             print(f"# bit {i}")
    #             p = cand.preimages(w)

    #             i_zeros = set(k for k, v in cand.map.items() if i in v)

    #             print(f"  {i_zeros = }")

    #             could_map_zero = []
    #             for x in p:
    #                 x_zeros = set(j for j in range(n) if (perm[x] & (1 << j)) == 0)
    #                 print(f"   {x=} {x_zeros=}")

    #                 z = x_zeros & i_zeros

    #                 if z and len(z) >= skipped:
    #                     could_map_zero.append(z)

    #             print(f"cmz = {could_map_zero} (len {len(could_map_zero)})")
    #             done = False
    #             c = None
    #             removed = 0
    #             if len(could_map_zero) == 0:
    #                 # print(f"!! Early remove: {zz} -//-> {w}")
    #                 # this is a "one" position. none of the zeros should map
    #                 # here
    #                 pass

    #             if len(could_map_zero) == 1:
    #                 done = True
    #                 c = could_map_zero[0]

    #             # if multiple matches, take the union
    #             # this is NOT the optimal solution: nominally possible to solve
    #             # for which bit assignments are consistent with 
    #             elif len(could_map_zero) > 1:
    #                 skipped += 1
    #                 comb = reduce(operator.or_, could_map_zero)
    #                 assert(len(comb) > 0)
    #                 c = comb
    #                 done = False

    #             if c is not None:
    #                 if any(len(cand.map[c_]) > 1 for c_ in c): 
    #                     removed += len(c)
                    
    #                     for j in range(n):
    #                         if j not in c:
    #                             cand.map[j].discard(i)
    #                 else:
    #                     done = False

    #                 # TODO: not sure about this part
    #                 if len(c) == 1:
    #                     c0 = next(iter(c))

    #                     if w.bit_count() == 1:
    #                         wb = (w - 1).bit_length()

    #                         y = set()
    #                         z = set()
    #                         for x in p:
    #                             if (perm[x] & (1 << c0)) == 0:
    #                                 y.add(x)
    #                             else:
    #                                 z.add(x)

    #                         if len(y) == 1:
    #                             x = y.pop()
    #                             xb = (x - 1).bit_length()
    #                             print(f"!! Enforce: {x} -> {w}")
    #                             cand.map[xb] = {wb}

    #                         for zz in z:
    #                             print(f"!! Remove: {zz} -//-> {w}")
    #                             zb = (zz - 1).bit_length()
    #                             cand.map[zb].discard(wb)
    #                 else:
    #                     done = False
                        

    #                 # for k, v in cand.map.items():
    #                 #     if i not in v:
    #                 #         # bit k cannot map to i
    #                 #         continue

    #                 #     if perm[x] & (1 << k) == 0:
    #                 #         c2[k].add(i)
    #                 #         has_zero_at_i += 1
    #             # print("@", i, x, has_zero_at_i)
    #             # if has_zero_at_i == 1:
    #             #     break

    #             print(cand)
    #             if cand.normalize():
    #                 print("Normalized", cand, sep='\n')

    #             if done:
    #                 break

    #         # print(c2)
    #         # print(f"{has_zero_at_i} could have zero at b{i}")

    #             # for k, v in cand.map.items():
    #             #     cand.map[k] -= c2[k]


    #             # # no mapping gaives a zero at set i.
    #             # # move on to the next bit
    #             # if len(has_zero_at_i) == 0:
    #             #     continue

    #             # for j in range(n):
    #             #     if i == j:
    #             #         cand.map[i] &= has_zero_at_i
    #             #     else:
    #             #         cand.map[j] -= has_zero_at_i

    #             # if len(c2) == 1:
    #             #     break

    #             # possible mapping now is w -> z, for all z in the has_zero set

    #         print("=== Cand Updated: ===")
    #         print(cand)

    #         if cand.complete():
    #             break
        
    #     if cand.complete():
    #         break
        
    # if not cand.complete():
    #     os.abort()
        
    # return [cand.map[i].pop() for i in range(n)]

            # # and find the lex min.
            # print(f" Check bit {i}")
            # q = set()

            # # which mappings would send a zero to bit i?
            

            # # intersect mapping that sends w to x

            # for x in p:
            #     # given x -> w, does bit i = 0?
            #     ##### 
            #     if perm[x] & (1 << i) == 0:
            #         q.add(x)
            #     # print(f"  {x} bit {i} map to 0?")
            # print(q)

            # if len(q) == 0:
            #     # no 0s in this bit. all shuffles will have 1.
            #     # move on to next bit (i -= 1)
            #     continue

            # at least one 1 in this bit

q = canonicalize(perm)
pp = apply_shuffle(perm, q)
print(perm, q, pp)
assert(set(q) == set(range(n)))
assert(set(pp) == set(range(N)))

good = 0

for _ in range(128):
    shuffle(q)
    p2 = apply_shuffle(perm, q)

    if np.array_equal(perm, p2):
        continue

    q2 = canonicalize(p2)
    pp2 = apply_shuffle(perm, q2)

    good += 1

    print(pp)
    print(pp2)
    assert(np.array_equal(pp, pp2))


print(f"{good=}")

# s0 = next(G)

# i0 = next(s0)
# p0 = perm[i0]
# w = p0.bit_count()
# adj = 2 ** w - 1
# print(f"i={i0} -> {p0}; b{p0:0{n}b} = {adj}")

# for i in range(n):
#     c = cand.a[i, :]
#     # if this bit was set
#     if p0 & (1 << i):
#         # cannot map to a high order bit
#         # must map to LSB
#         c[w:] = False
#     else:
#         # zero maps to MSB
#         c[:w] = False

# print(cand)

# sz = next(G)
# iz = next(sz)
# pz = perm[iz]
# print(f"i={iz} -> {pz}; b{pz:0{n}b}")

# # Try to get a zero in the MSB
# msb_possible = np.where(cand.a[:, -1])[0]
# print(f"{msb_possible=}")
# msb_bits = np.array([(pz & 1 << m) != 0 for m in msb_possible])
# print("msb bits", np.array(msb_bits, dtype=int))

# # num_zeros = msb_bits.count(0)
# could_be_msb = msb_possible[msb_bits == False]
# couldnt_be_msb = msb_possible[msb_bits == True]
# print(f"{could_be_msb=}")
# print(f"{couldnt_be_msb=}")

# if len(could_be_msb):
#     cand.a[could_be_msb, :-len(could_be_msb)] = False
#     cand.a[couldnt_be_msb, -len(could_be_msb):] = False

# print(cand)


# print("== Checking LSB chunk ==")
# # no selectivity from MSB chunk. Consider LSB chunk
# # technically also do this if both all MSB chunk zero
# lsb_possible = np.where(cand.a[:, 0])[0]
# print(f"{lsb_possible=}")
# lsb_bits = np.array([(pz & 1 << m) != 0 for m in lsb_possible])

# print("lsb chunk bits", np.array(lsb_bits, dtype=int))
# # now look at *msb* of lsb chunk
# could_be_msb = lsb_possible[lsb_bits == False]
# couldnt_be_msb = lsb_possible[lsb_bits == True]

# print(f"{could_be_msb=}")
# print(f"{couldnt_be_msb=}")

# if len(could_be_msb):
#     s = w
#     cand.a[could_be_msb, 0:s-len(could_be_msb):] = False
#     cand.a[couldnt_be_msb, s-len(could_be_msb):s] = False


# ## vaguely
# # divide into chunks (unique partitions)
# # for each chunk, in order (L to R, MSB to LSB), consider input bits at those
# # positions. any zeros? enforce that swap. if no swap happened, check next
# # chunk (this is more lexicographical)


# print(cand)

# s2 = list(next(G))

# for k in s2:
#     v = perm[k]
#     print(f"{k:4} -> {v:4} b{v:0{n}b}")


# # could_become_p1 = np.where(cand.a[:, 0])[0]
# # print(f"{could_become_p1}")

# # for each bit, find optimal assignment for perm(2^index)
# for index in range(n):
#     index_pre = np.where(cand.a[:, index])[0]
#     print(f"///// perm {index=} -> π'({1 << index}) could come from π({1 << index_pre})")
#     print("Start at MSB.")
#     for b in range(n - 1, 0 - 1, -1):
#         index_pre = np.where(cand.a[:, index])[0]
#         # where could this bit come from, given the current candidate mapping?
#         preimage = np.where(cand.a[:, b])[0]

#         print("output bit", b, "has preimage idx", preimage)

#         # which powers of two have zeros in one of these locations?
#         zeros = set()
#         zero_idx = set()
#         for i in index_pre:
#             p = perm[1 << i]
#             r = set(int(q) for q in preimage if p & (1 << q) == 0)
#             if r:
#                 zeros.add(i)
#                 zero_idx.update(r)
#             else:
#                 print(f"  perm idx={i} contains no zero in bit preimage")
#                 # cand.a[b, i] = False

#         print("valid perms with zero in any position:", 1 << np.array(list(zeros)) if zeros else "{}")
#         print(f"those bit positions: {zero_idx}")

#         if len(zeros) == 1 and ((zeros == zero_idx and b != index) or (zeros != zero_idx and index == b)):
#             print("Minimal not possible - skip")
#             continue

#         marked = False

#         if len(zeros) == 1:
#             z = zeros.pop()
#             if sum(cand.a[z, :]) > 1:
#                 print(f"Marking (1) {z} -> {index}")
#                 marked = True
#                 cand.a[:, index] = False
#                 cand.a[z, :] = False
#                 cand.a[z, index] = True
#         elif 1 < len(zeros) < len(index_pre):
#             print("==> FILTER z")
#             for i in index_pre:
#                 if i not in zeros:
#                     if sum(cand.a[i, :]) > 1:
#                         cand.a[i, index] = False

        
#         if not marked:
#             if len(zero_idx) == 1:
#                 z = zero_idx.pop()
#                 print(f"Marking (2) {z} -> {b}")
#                 cand.a[:, b] = False
#                 cand.a[z, :] = False
#                 cand.a[z, b] = True
#             elif 0 < len(zero_idx) < len(preimage):
#                 print("==> FILTER zi")
                
#                 for i in preimage:
#                     if i not in zero_idx:
#                         if sum(cand.a[i, :]) > 1:
#                             cand.a[i, b] = False

#                         # col_sum = np.sum(cand.a, axis=1)
#                         # row_sum = np.sum(cand.a, axis=0)

#                         # if any(col_sum == 0) or any(row_sum) == 0:
#                         #     # revert!
#                         #     print("oops, revert")
#                         #     cand.a[i, b] = True
#                     else:
#                         pass
#                         # if (sum(cand.a[i, :]) > 1):
#                         #     cand.a[i, :] = False
#                         #     cand.a[i, b] = True

#         print(cand)

#         old = cand.a

#         if cand.complete():
#             break

#         cand.normalize()
#         # print("Normalized")
#         # print(cand)


#     if cand.complete():
#         print("DONE")
#         break

#     if not cand.legal():
#         print("ERROR!")
#         sys.exit(1)
#         break

#     print(cand)

# print("FINAL:")
# print(cand)

# groups = defaultdict(set)

# for i in range(n):
#     groups[tuple(cand.a[i])].add(i)

# print(f"{groups=}")

# for k, v in groups.items():
#     if k[-1] == True:
#         print("MSB group", np.array(k, dtype=int))

#     else:
#         print("LSB group...", k)



# l_group = cand.a[0]
# for i in range(n):
#     r_group = cand.a[i]
#     if not all(r_group == l_group):
#         break

# print("L", l_group, "R", r_group)

# min_l = np.dot((2 ** sum(l_group & pz) - 1), bincomp)
# min_r = np.dot((2 ** sum(r_group & pz) - 1), bincomp)

# print(min_l, min_r)

# mapped = (min_l) << w | min_r 
# print("mapped:", mapped)

# msb_chunk_from = np.where(cand.a[:, n - 1])[0]
# msb_chunk_size = n - w
# print(msb_chunk_from, msb_chunk_size)

# q = reduce(operator.or_, 1 << msb_chunk_from)
# print(bin(q & pz))

# msb_weight = (q & pz).bit_count()
# # shift these over
# print(msb_weight)

# if msb_weight < msb_chunk_size:
#     # do some filtering
#     for i, m in enumerate(msb_chunk_from):
#         if i < msb_weight:
#             s = n - msb_chunk_size
#             cand[:, s:(s+msb_weight)] = False
#         else:
#             s = n - msb_weight
#             cand[:, -s] = False
#         print(f"{m=}")


#     print(f"{msb_weight=}")
#     print(f"{msb_chunk_size=}")
# else:
#     print("MSB chunk all ones, no filter")

# print(cand)

# # for each bit in p0: remove zero locations
# # for each zero in p0: remove one locations