import matplotlib.pyplot as plt
import circuit_lib as lib
import random

seed = random.randint(0, 1e6)
random.seed(seed)
print(seed)

# c = lib.random_circuit(5, 4)
c = lib.Circuit([
    lib.Gate(3, 2, 5),
    lib.Gate(3, 1, 5),
    lib.Gate(0, 2, 5),
    lib.Gate(0, 1, 4)])
# c = lib.Circuit([
#     lib.Gate(0, 2),
#     lib.Gate(0, 2),
#     lib.Gate(1, 2),
#     lib.Gate(1, 3)
# ])
print(c)

g = lib.SkeletonGraph.from_circuit(c, True)

# g.draw()

# for n in g.nodes:
#     print(f"=== Node {n}")
#     print(" in ", [e['weight'] for (_, _, e) in g.in_edges(n, data=True)])
#     print(" out", [e['weight'] for (_, _, e) in g.out_edges(n, data=True)])
# plt.show()

g.canonical()

# plt.show()