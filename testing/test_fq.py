import pytest
from field.fq import Fq


def test_fq_without_field_modulus():
    with pytest.raises(AttributeError):
        a = Fq.zero()


def test_fq_with_invalid_val():
    Fq.field_modulus = 17
    with pytest.raises(TypeError):
        a = Fq("Yolo")


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
    c = a + b
    assert(c == b)
    assert(c.val == 1)


def test_fq_sub():
    Fq.field_modulus = 17
    a = Fq.one()
    b = Fq.one()
    res = Fq.zero()
    c = a - b
    assert(c == res)


def test_fq_mul():
    Fq.field_modulus = 17
    a = Fq(3)
    b = Fq(2)
    res = Fq(6)
    c = a * b
    assert(c == res)


def test_fq_neg():
    Fq.field_modulus = 17
    a = -Fq(3)
    res = Fq(14)
    assert(a == res)


def test_fq_div():
    Fq.field_modulus = 17
    a = Fq(2)
    b = Fq(9)
    res = Fq(4)
    c = a / b
    assert(c == res)


def test_fq_exp():
    Fq.field_modulus = 17
    a = Fq(3) ** 2
    res = Fq(9)
    assert(a == res)