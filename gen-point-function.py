#!/usr/bin/env python3
#
# Generate a point function circuit with an n-bit key, over (n+2) wires.
import argparse
import sys
import math

parser = argparse.ArgumentParser(description="Generate point function circuit with an n-bit key.")
parser.add_argument('-n', type=int, required=True, help='num bits')
parser.add_argument('-k', type=int, required=True, help='secret key')

args = parser.parse_args()

n = args.n
max_k = (1 << n) - 1
if args.k < 0 or args.k > max_k:
    print(f"Error: k must be between 0 and {max_k} (inclusive) for {n=}.")
    sys.exit(1)

bits = [i for i in range(n) if args.k & (1 << i) == 0]

# flip zeros into ones
for b in bits:
    print(f"! {b}")

def generate_tof(active, controls):
    empty = list(set(range(n + 2)) - set(controls) - {active})
    # print(f"Gen Tof[{active}] ctrls {controls}")
    # print("  ", empty)

    base_gate = f"{active} {controls[0]} {empty[0]}"

    stair = []
    last_tof = None
    for i in range(1, len(controls) - 2):
        stair.append(f"{empty[i - 1]} {controls[i]} {empty[i]}")
        last_tof = empty[i]
    
    print(base_gate)
    print('\n'.join(stair))
    print(f"{last_tof} {controls[-2]} {controls[-1]}")
    print('\n'.join(reversed(stair)))
    print(base_gate)
    print('\n'.join(stair))
    print(f"{last_tof} {controls[-2]} {controls[-1]}")
    print('\n'.join(reversed(stair)))


# first decomposed Tof has active=n+2 and ctrl [0, n/2]
# second decomposed Tof has active=n+2 and ctrl [n/2, n]
# ancilla n+1

b = math.ceil(n/2)

generate_tof(n, list(range(0, b)))
generate_tof(n+1, list(range(b, n + 1)))
generate_tof(n, list(range(0, b)))
generate_tof(n+1, list(range(b, n + 1)))

# flip back
for b in bits:
    print(f"! {b}")

