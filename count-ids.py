#!/usr/bin/env python3

import circuit_lib as lib
from rich.progress import Progress
from collections import defaultdict
import multiprocessing as mp
import os
import random
import pickle
from pathlib import Path
import gzip

n = 5
m = 8

def invert(perm):
    n = len(perm)
    inv = [0] * n

    for i, p in enumerate(perm):
        inv[p] = i

    return tuple(inv)

def ckt_work(gates):
    ckt = lib.Circuit(gates, wires=n)
    g = lib.SkeletonGraph.from_circuit(ckt)
    p = tuple(ckt.perm())
    s = tuple(g.sorted_gates())

    return (p, s)

def generate_new(n, m_):
    all_ckt = lib.all_circuits(n, m_)
    count = (n * (n * n - n)) ** m_

    perms = defaultdict(set)
    revmap = {}

    print(f"Checking {count} circuits...")

    with Progress() as prog:
        task = prog.add_task("Build circuits...", total=count)

        with mp.Pool(os.process_cpu_count() // 2) as pool:
            par = pool.imap_unordered(
                ckt_work, all_ckt
            )

            for c in par:
                p, s = c

                invp = invert(p)
        
                # If inverse already listed, then add this circuit there
                # First argument of tuple says whether this circuit is inverted or not
                if invp in perms:
                    perms[invp].add((False, s))
                    revmap[s] = (False, invp)
                else:
                    # otherwise, add to own perm
                    perms[p].add((True, s))
                    revmap[s] = (True, p)
                
                prog.advance(task)

    return perms, revmap

def main():
    # generate all circuits of this size
    assert (m % 2) == 0
    m_ = m // 2
    count = (n * (n * n - n)) ** m_

    picklename = f'../pickles/n{n}m{m_}.gz'
    if Path(picklename).exists():
        print("Load pickle...")
        perms, revmap = pickle.load(gzip.open(picklename, 'rb'))
    else:
        perms, revmap = generate_new(n, m_)
        pickle.dump((perms, revmap), gzip.open(picklename, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    singletons = 0

    maxp = None
    maxc = 0

    for p, v in perms.items():
        # print(f"{str(p):40} {len(v):2}")
        if len(v) == 1:
            singletons += 1

        if len(v) > maxc:
            maxp = p
            maxc = len(v)

        inv = invert(p)
        assert (p == inv) or (inv not in perms)

    print(len(perms), "perms")
    pct = len(perms) / count
    print(f"{count} circuits ({pct*100:.1f}%)")

    print(singletons, "singletons")
    print(f"Max perm ({maxc}): {maxp}")

    print("\n")

    rand_ckt = random.choice(list(revmap.keys()))
    print(lib.Circuit(rand_ckt))
    do_invert, p_ = revmap[rand_ckt]
    print("[Inverse]\n" if do_invert else "", end='')

    invp = invert(p_)
    if invp == p:
        print("Own inverse!")

    p = invp if do_invert else p_

    others = perms[p]
    print(len(others), "others:")
    for i, k in others:
        print("----")
        print("[Inverse]\n" if i else "", end='')
        print(lib.Circuit(k))




if __name__ == '__main__':
    main()