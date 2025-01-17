#!/usr/bin/env python3

from circuit_lib import SkeletonCache
import argparse
import itertools as it

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-n", type=int, help="number of wires/bits")
    parse.add_argument("-m", type=int, help="number of gates")
    parse.add_argument("-r", action="store_true", help="run all")
    args = parse.parse_args()


    if args.n and args.m:
        cache = SkeletonCache(args.n, args.m)

        cache.build()
        
        # cache.print()
        cache.stats()

    elif args.r:
        max_combined = 12

        skel_list = []

        for (n, m) in it.product(range(2, max_combined), range(2, max_combined)):
            if n + m > max_combined:
                continue
            cache = SkeletonCache(n, m)
            skel_list.append(cache)

        # sort by num circuits
        for cache in sorted(skel_list, key=lambda c: c.num_ckts):
            print("—" * 32)
            cache.build()
            cache.stats()

if __name__ == '__main__':
    main()