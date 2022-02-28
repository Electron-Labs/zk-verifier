import pytest
from field.fq2 import Fq2, Fq


def test_fq2_without_non_residue():
    with pytest.raises(AttributeError):
        a = Fq2(2)


def test_fq2_without_list():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    with pytest.raises(AssertionError):
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
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = Fq2([Fq(9), Fq(7)])
    b = Fq2([Fq(14), Fq(15)])
    res = Fq2([Fq(6), Fq(5)])
    c = a + b
    assert(c == res)


def test_fq2_sub():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = Fq2([Fq(9), Fq(7)])
    b = Fq2([Fq(14), Fq(15)])
    res = Fq2([Fq(12), Fq(9)])
    c = a - b
    assert(c == res)


def test_fq2_mul():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = Fq2([Fq(9), Fq(7)])
    b = Fq2([Fq(14), Fq(15)])
    res = Fq2([Fq(4), Fq(3)])
    c = a * b
    assert(c == res)


def test_fq2_neg():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = -Fq2([Fq(9), Fq(7)])
    res = Fq2([Fq(8), Fq(10)])
    assert(a == res)


def test_fq2_inverse():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    a = Fq2([Fq(9), Fq(7)])
    b = a.inverse()
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