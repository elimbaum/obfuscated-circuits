import random
import circuit_lib as lib
import matplotlib.pyplot as plt
from pprint import pprint
from collections import Counter, defaultdict
from functools import partial
from rich.progress import track
from rich import print
from rich.pretty import pprint


def random_circuit(n, m):
    b = list(lib.base_perms(n))
    return lib.Circuit(random.choices(b, k=m))

ckt = random_circuit(6, 6)

print(ckt)
g = lib.SkeletonGraph.from_circuit(ckt, augmented=True)
print(g)

# g.draw()

order = 0

dp = None

dp = g.deg_profile()

# plt.show()

input()

def ckt_to_prof(n, c):
    ckt = lib.Circuit(c, wires=n)
    g = lib.SkeletonGraph.from_circuit(ckt, augmented=True)
    return g.canonical()

n = 4
m = 4
dc = defaultdict(list)

print(f"all: {n=}, {m=}")

b = len(list(lib.base_perms(n)))

total = b ** (m - 1)
total *= len(list(lib.zero_active_perms(n)))
n_ckt = b ** m

for c in track(lib.all_circuits(n, m, True), total=total):
    dc[ckt_to_prof(n, c)].append(c)
# dc.update({ckt_to_prof(n, c): c for c in lib.all_circuits(n, m)})
# pprint(dc)

for sk, cs in dc.items():
    # print("\n\n==== Skeleton ====")
    # print(sk)
    # print(cs)
    gc = None
    for c in cs:

        gc_ = Counter(map(lambda g: g.type.value, c))
        if gc is None:
            gc = gc_
        elif gc != gc_:
            print("ERROR!")
            print("  ", end='')
            print(c)
            print(cs)
            print(gc)
            print(gc_)
            # Structurally equivalent circuits should have the same gates.
            assert(gc == gc_)

print(f"{n=} {m=} {len(dc)=} {len(dc) / n_ckt * 100:.1f}%")