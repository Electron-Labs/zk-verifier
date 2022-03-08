import pytest
from bn128.bn128 import BN128
from bn128.g1 import G1
from bn128.g2 import G2
from field.fq12 import Fq12

def test_bn128_pairing():
    bn128 = BN128()
    a = 40
    b = 75

    g1a = bn128.g1.mul_scalar(a)
    g2a = bn128.g2.mul_scalar(b)

    g1b = bn128.g1.mul_scalar(b)
    g2b = bn128.g2.mul_scalar(a)

    pre1a = bn128.ate_precompute_g1(g1a)
    pre2a = bn128.ate_precompute_g2(g2a)

    pre1b = bn128.ate_precompute_g1(g1b)
    pre2b = bn128.ate_precompute_g2(g2b)

    r1 = bn128.ate_miller_loop(pre1a, pre2a)
    r2 = bn128.ate_miller_loop(pre1b, pre2b)

    rbe = r1 * r2.inverse()
    res = bn128.final_exponentation(rbe)

    assert(isinstance(res, Fq12))

    one = Fq12.one()

    assert(one == res)