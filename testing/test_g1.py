import pytest
from field.fq import Fq
from bn128.g1 import G1


def test_g1_add():
    G1.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583

    gr1 = G1([1, 2])
    gr2 = G1([1, 2])

    r1 = Fq(33)
    r2 = Fq(44)

    gr1 = gr1.mul_scalar(r1)
    gr2 = gr2.mul_scalar(r2)

    r1r2 = r1 + r2

    grsum1 = gr1 + gr2
    grsum2 = G1([1, 2]).mul_scalar(r1r2)

    assert(grsum1 == grsum2)

    a = grsum1.affine()

    assert(a.val[0] == 21526464323725832663882905544083280657770325585151797133383551854196089356032)
    assert(a.val[1] == 8545759555567142326482563981456384114560528812235279370284511019768507753138)


def test_g1_double():
    G1.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583

    gr1 = G1([1, 2]).mul_scalar(2)
    gr2 = G1([1, 2]).double()

    assert(gr1 == gr2)