# test permutations

import perm_canon2 as pc
import random
import numpy as np

seed = 0
random.seed(0)
print(f"{seed=}")

REPS = 16

SCALE_FACTOR = 5

def consistent_canon(perm, n):
    N = 1 << n
    q = pc.canonicalize(perm, do_print=True)
    assert set(q) == set(range(n))
    
    k = pc.apply_shuffle(perm, q)
    assert set(k) == set(range(N))

    for _ in range(REPS):
        # get a new bit shuf
        random.shuffle(q)

        pp = pc.apply_shuffle(perm, q)

        q2 = pc.canonicalize(pp, do_print=True)
        k2 = pc.apply_shuffle(pp, q2)

        assert np.array_equal(k, k2)

def random_helper(n):
    N = 1 << n
    perm = np.arange(N)
    random.shuffle(perm)

    consistent_canon(perm, n)

def test_three():
    for _ in range(SCALE_FACTOR * 10):
        random_helper(3)

def test_four():
    for _ in range(SCALE_FACTOR * 10):
        random_helper(4)

def test_five():
    for _ in range(SCALE_FACTOR * 10):
        random_helper(5)

def test_six():
    for _ in range(SCALE_FACTOR * 5):
        random_helper(6)

def test_ten():
    for _ in range(SCALE_FACTOR):
        random_helper(10)

def test_ten():
    for _ in range(SCALE_FACTOR):
        random_helper(12)

def test_identity():
    n = 6
    idperm = np.arange(1 << n)
    q = pc.canonicalize(idperm, False)
    assert np.array_equal(q, np.arange(n))

def test_near_identity():
    n = 6
    p = np.arange(1 << n)
    s = list(pc.strings_of_weight(n // 2, n))

    # swap some heavy entries
    p[s[-1]], p[s[-2]] = p[s[-2]], p[s[-1]]

    print(p)

    consistent_canon(p, n)