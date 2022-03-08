import pytest
from bn128.bn128 import BN128
from bn128.g1 import G1
from bn128.g2 import G2

def test_bn128_pairing():
    bn128 = BN128()