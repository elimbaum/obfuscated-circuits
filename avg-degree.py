import random
import circuit_lib as lib
import matplotlib.pyplot as plt
from pprint import pprint
from collections import Counter, defaultdict
from functools import partial
from rich.progress import track
from rich import print
from rich.pretty import pprint
import numpy as np

def avg_graph_weight(n, m, gates, trials=2000):
    
    def random_tof_ckt():
        return lib.Circuit(random.choices(gates, k=m))
    
    # expected = 2 * m * m / n
    expected = m*(m-1)/2 * (4*(n-1))/(n*n)

    dd = []
    for i in track(range(trials)):
        ckt = random_tof_ckt()
        g = lib.SkeletonGraph.from_circuit(ckt, augmented=False)
        dd.append(g.number_of_edges())

    actual = np.mean(dd)
    pct = actual / expected
    print(f"{n=},{m=}: {expected=}, got {actual:.2f} ({pct*100:.1f}%)")
    if actual > expected:
        print("!!! above bound")

# tall circuits
for n in [50, 100, 200, 500]:
    print(f"----- {n=}")
    gates = list(filter(lambda g: g.type == lib.GateT.CCNOT, lib.base_perms(n)))
    for m in [2, 5, 10, 20, 50, 100, 200, 500]:
        avg_graph_weight(n, m, gates)