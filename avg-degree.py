# %%
import random
# import circuit_lib as lib
import matplotlib.pyplot as plt
from pprint import pprint
from collections import Counter, defaultdict
from functools import partial
from rich.progress import track
from rich import print
from rich.pretty import pprint
import numpy as np
import itertools as it

# %%

# def avg_graph_weight(n, m, gates, trials=2000):
    
#     def random_tof_ckt():
#         return lib.Circuit(random.choices(gates, k=m))
    
#     # expected = 2 * m * m / n
#     expected = m*(m-1)/2 * (4*(n-1))/(n*n)

#     dd = []
#     for i in track(range(trials)):
#         ckt = random_tof_ckt()
#         g = lib.SkeletonGraph.from_circuit(ckt, augmented=False)
#         dd.append(g.number_of_edges())

#     actual = np.mean(dd)
#     pct = actual / expected
#     print(f"{n=},{m=}: {expected=}, got {actual:.2f} ({pct*100:.1f}%)")
#     if actual > expected:
#         print("!!! above bound")

# # tall circuits
# for n in [50, 100, 200, 500]:
#     print(f"----- {n=}")
#     gates = list(filter(lambda g: g.type == lib.GateT.CCNOT, lib.base_perms(n)))
#     for m in [2, 5, 10, 20, 50, 100, 200, 500]:
#         avg_graph_weight(n, m, gates)

n = 4
m = 50

# [active, ctrl, ctrl]
gates = np.array([np.random.choice(np.arange(n), 3, replace=False) for _ in range(m)])
print(gates)

# %%
def partition_runs(m):
    m = np.array(m)
    nz = m[m.nonzero()]
    if not len(nz):
        return []

    sign = np.sign(nz)
    split = np.flatnonzero(sign[1:] != sign[:-1]) + 1
    return [np.abs(x).tolist() for x in np.split(nz, split)]

# %%
active = np.equal.outer(gates[:, 0], np.arange(n))
print(active)

ctrl1 = np.equal.outer(gates[:, 1], np.arange(n))
ctrl2 = np.equal.outer(gates[:, 2], np.arange(n))
print(ctrl1 | ctrl2)
print(active - (ctrl1 | ctrl2).astype(int))

z = active - (ctrl1 | ctrl2).astype(int)
# %%
q = z * np.c_[np.arange(1, m + 1)]
q

# %%
collisions = set()

for c in q.T:
    p = partition_runs(c)
    print(p)
    for h in it.pairwise(p):
        collisions.update(frozenset(q) for q in it.product(*h))
        
# %%
collisions
# %%
len(collisions)
# %%


def simulate(n, m):
    ckt = np.array([np.random.choice(np.arange(n), 3, replace=False) for _ in range(m)])
    active = np.equal.outer(ckt[:, 0], np.arange(n))

    ctrl1 = np.equal.outer(ckt[:, 1], np.arange(n))
    ctrl2 = np.equal.outer(ckt[:, 2], np.arange(n))

    z = active - (ctrl1 | ctrl2).astype(int)
    q = z * np.c_[np.arange(1, m+1)]

    collisions = set()
    
    for c in q.T:
        p = partition_runs(c)
        for h in it.pairwise(p):
            collisions.update(frozenset(q) for q in it.product(*h))
    
    return len(collisions)

# %%
for N in [4, 6, 8, 16, 32]:
# N = 4
    print(f"{N=}")
    mm = []
    L = []
    for t in range(10):
        # print(t)
        for m in range(40, 2000, 40):
            m += np.random.randint(-20, 20)
            mm.append(m)
            L.append(simulate(N, m))

    # density = L / (mm * (np.array(mm) - 1) / 2)
    density = 2 * np.array(L) / mm

    P = plt.scatter(mm, density, alpha=0.5, s=4)
    # x = np.arange(10, 1000)
    sc = 6 + 6 * (1 - 2/N)**3
    sc1 = 6 + (6 - 9/(3 * N/2))
    sc2 = 12 * (1 - 1 / N)
    # plt.hlines(sc, 10, 2000, linestyles="--", color=P.get_facecolor()[0], alpha=1)
    # plt.hlines(sc1, 10, 2000, linestyles="-.", color=P.get_facecolor()[0], alpha=1)
    plt.hlines(sc2, 10, 2000, linestyles=":", color=P.get_facecolor()[0], alpha=1)
    # plt.plot(x, sc, color=P.get_facecolor()[0])
    # plt.ylim(0, 0.1)

x = np.arange(10, 1000)
plt.hlines(12, 10, 2000, color='black')

plt.ylim(7, 13)
# %%
plt.show()
# %%
density
# %%
