#!/usr/bin/env python3

from circuit_lib import SkeletonCache
import argparse

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


        last_n = n = m = 2
        while n + m < 12:
            print("—" * 32)
            cache = SkeletonCache(n, m)
            cache.build()
            cache.stats()
            n -= 1
            m += 1

            if n < 2:
                last_n += 1
                n = last_n
                m = 2


if __name__ == '__main__':
    main()