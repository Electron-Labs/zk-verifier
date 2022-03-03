import pytest
from field.fq2 import Fq2
from field.fq import Fq


def test_fq2_without_non_residue():
    Fq2.non_residue = None
    with pytest.raises(AttributeError):
        a = Fq2(2)


def test_fq2_with_smaller_list():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    with pytest.raises(AssertionError):
        a = Fq2([2])


def test_fq2_with_bigger_list():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    with pytest.raises(AssertionError):
        a = Fq2([2, 3, 4])


def test_fq2_zero():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    zero = Fq2.zero()
    for i in range(Fq2.degree):
        assert(zero.val[i] == Fq.zero())


def test_fq2_one():
    Fq2.non_residue = Fq(16)
    Fq.field_modulus = 17
    one = Fq2.one()
    assert(one.val[0] == Fq.one())
    for i in range(1, Fq2.degree):
        assert(one.val[i] == Fq.zero())


def test_fq2_add():
    Fq.field_modulus = 7
    Fq2.non_residue = Fq(-1)
    a = Fq2([Fq(4), Fq(4)])
    b = Fq2([Fq(3), Fq(4)])
    res = Fq2([Fq(0), Fq(1)])
    c = a + b
    assert(c == res)


def test_fq2_sub():
    Fq.field_modulus = 7
    Fq2.non_residue = Fq(-1)
    a = Fq2([Fq(5), Fq(3)])
    b = Fq2([Fq(7), Fq(2)])
    res = Fq2([Fq(5), Fq(1)])
    c = a - b
    assert(c == res)


def test_fq2_mul():
    Fq.field_modulus = 7
    Fq2.non_residue = Fq(-1)
    a = Fq2([Fq(4), Fq(4)])
    b = Fq2([Fq(3), Fq(4)])
    res = Fq2([Fq(3), Fq(0)])
    c = a * b
    assert(c == res)


def test_fq2_neg():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = -Fq2([Fq(9), Fq(7)])
    res = Fq2([Fq(8), Fq(10)])
    assert(a == res)


def test_fq2_square():
    Fq.field_modulus = 7
    Fq2.non_residue = Fq(-1)

    a = Fq2([Fq(4), Fq(4)])
    b = a.square()
    c = a * a
    assert(b == c)

def test_fq2_inverse():
    Fq.field_modulus = 7
    Fq2.non_residue = Fq(-1)
    a = Fq2([Fq(4), Fq(4)])
    c = Fq2([Fq(1), Fq(6)])
    b = a.inverse()
    assert(b == c)
    # Verify basic property of inverse a*a^-1 = 1 mod p
    c = a * b
    res = Fq2.one()
    assert(c == res)


def test_fq2_mul_scalar():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = Fq2([Fq(2), Fq(3)])
    b = Fq(3)
    c = a.mul_scalar(b)
    res = Fq2([Fq(6), Fq(9)])
    assert(c == res)


def test_fq2_frobenius_map():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    a = Fq2([Fq(2), Fq(3)])
    a.frobenius_coeffs_c1 = [Fq(1),
                             Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)]
    a = a.frobenius_map(3)
    res = Fq2([Fq(2), Fq(21888242871839275222246405745257275088696311157297823662689037894645226208580)])
    assert(a == res)