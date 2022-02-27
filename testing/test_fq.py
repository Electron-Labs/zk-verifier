import pytest
from field.fq import Fq


def test_fq_without_field_modulus():
    with pytest.raises(Exception):
        a = Fq.zero()


def test_fq_zero():
    Fq.field_modulus = 17
    a = Fq.zero()
    assert(a.val == 0)
    assert(isinstance(a, Fq))


def test_fq_one():
    Fq.field_modulus = 17
    a = Fq.one()
    assert(a.val == 1)
    assert(isinstance(a, Fq))


def test_fq_add():
    Fq.field_modulus = 17
    a = Fq.zero()
    b = Fq.one()
    a.add(b)
    assert(a.eq(b))
    assert(a.val == 1)


def test_fq_sub():
    Fq.field_modulus = 17
    a = Fq.one()
    b = Fq.one()
    c = Fq.zero()
    a.sub(b)
    assert(c.eq(a))


def test_fq_mul():
    Fq.field_modulus = 17
    a = Fq(3)
    b = Fq(2)
    c = Fq(6)
    a.mul(b)
    assert(a.eq(c))


def test_fq_mul_scalar():
    Fq.field_modulus = 17
    a = Fq(3)
    b = 2
    c = Fq(6)
    a.mul_scalar(b)
    assert(a.eq(c))


def test_fq_neg():
    Fq.field_modulus = 17
    a = Fq(3)
    b = Fq(14)
    a.neg()
    assert(a.eq(b))


def test_fq_inverse():
    Fq.field_modulus = 17
    a = Fq(2)
    b = Fq(9)
    a.inverse()
    assert(a.eq(b))


def test_fq_div():
    Fq.field_modulus = 17
    a = Fq(2)
    b = Fq(9)
    c = Fq(4)
    a.div(b)
    assert(c.eq(a))
    assert(b.eq(b))


def test_fq_square():
    Fq.field_modulus = 17
    a = Fq(3)
    b = Fq(9)
    a.square()
    assert(a.eq(b))
