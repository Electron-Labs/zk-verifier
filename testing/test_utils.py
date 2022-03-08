from field.utils import find_naf
import random


def is_naf_is_valid(num, naf):
    base = 0
    idx = 0
    for n in naf:
        if n == -1:
            base -= (2 ** idx)
        elif n == 1:
            base += (2 ** idx)
        idx += 1
    return num == base


def test_naf_representation():
    randomlist = random.sample(range(100, 42355), 1000)
    for r in randomlist:
        naf = find_naf(r)
        assert(is_naf_is_valid(r, naf))