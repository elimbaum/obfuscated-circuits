#!/usr/bin/env python3

from circuit_lib import SkeletonCache
import argparse

def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("-n", type=int, help="number of wires/bits")
    parse.add_argument("-m", type=int, help="number of gates")
    args = parse.parse_args()

    cache = SkeletonCache(args.n, args.m)

    cache.build()
    
    # cache.print()
    cache.stats()

if __name__ == '__main__':
    main()