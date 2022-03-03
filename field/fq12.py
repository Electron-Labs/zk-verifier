from field.fq import Fq
from field.fq2 import Fq2
from field.fq6 import Fq6

class Fq12:
    non_residue = None
    degree = 12
    frobenius_coeffs_c1 = None

    def __init__(self, val):
        if not self.non_residue:
            raise AttributeError("Quadratic non-residue hasn't been specified")

        assert(isinstance(self.non_residue, Fq2))

        self.val = val

    @classmethod
    def zero(cls):
        return cls([Fq6.zero(), Fq6.zero()])

    @classmethod
    def one(cls):
        return cls([Fq6.one(), Fq6.zero()])

    def __add__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x + y for x, y in zip(self.val, other.val)])

    def __sub__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x - y for x, y in zip(self.val, other.val)])

    def __neg__(self):
        return type(self)([-c for c in self.val])

    def __eq__(self, other):
        assert isinstance(other, self.__class__)
        for c1, c2 in zip(self.val, other.val):
            if c1 != c2:
                return False
        return True

    def mul_by_residue(self, val):
        return Fq6([val.val[2] * self.non_residue, val.val[0], val.val[1]])

    def __mul__(self, other):
        # Multiplication and Squaring on Pairing-Friendly.pdf; Section 3 (Karatsuba)
        # https://pdfs.semanticscholar.org/3e01/de88d7428076b2547b60072088507d881bf1.pdf

        a0 = self.val[0]
        a1 = self.val[1]
        b0 = other.val[0]
        b1 = other.val[1]

        v0 = a0 * b0
        v1 = a1 * b1

        val1 = v0 + self.mul_by_residue(v1)
        val2 = (a0 + a1) * (b0 + b1) - (v0 + v1)

        return self.__class__([val1, val2])

    def inverse(self):
        # High-Speed Software Implementation of the Optimal Ate Pairing over Barreto–Naehrig Curves .pdf
        # https://eprint.iacr.org/2010/354.pdf , algorithm 8
        t0 = self.val[0] * self.val[0]
        t1 = self.val[1] * self.val[1]
        t2 = t0 - self.mul_by_residue(t1)
        t3 = t2.inverse()
        return self.__class__([(self.val[0] * t3), -(self.val[1] * t3)])

    def __div__(self, other):
        return self * other.inverse()

    def __truediv__(self, other):
        return self.__div__(other)

    def __repr__(self):
        return repr(self.val)

    def mul_scalar(self, base):
        res = Fq12.zero()
        rem = Fq(base)
        exp = self

        while rem.val != 0:
            if rem.val % 2 == 1:
                res = res + exp

            exp += exp
            rem.val= rem.val >> 1

        return res

    def __pow__(self, other):
        if other == 0:
            return Fq12.one()
        elif other == 1:
            return Fq12(self.val)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    def frobenius_map(self, power):
        assert(self.frobenius_coeffs_c1 is not None)
        return self.__class__([
            self.val[0].frobenius_map(power),
            self.val[1].frobenius_map(power).mul_by_fq2(self.frobenius_coeffs_c1[power % 12]),
        ])