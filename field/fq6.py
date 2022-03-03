from field.fq import Fq
from field.fq2 import Fq2

class Fq6:
    non_residue = None
    degree = 6
    frobenius_coeffs_c1 = None
    frobenius_coeffs_c2 = None

    def __init__(self, val):
        if not self.non_residue:
            raise AttributeError("Quadratic non-residue hasn't been specified")

        assert(isinstance(self.non_residue, Fq2))

        self.val = val

        assert(len(self.val) > 0)
        inner_len = len(self.val[0].val)
        assert(inner_len == 2)
        assert(isinstance(self.val[0], Fq2))
        assert(len(self.val) == self.degree // inner_len)

    @classmethod
    def zero(cls):
        return cls([Fq2.zero(), Fq2.zero(), Fq2.zero()])

    @classmethod
    def one(cls):
        return cls([Fq2.one(), Fq2.zero(), Fq2.zero()])

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
    
    def __mul__(self, other):
        a0, b0 = self.val[0], other.val[0]
        a1, b1 = self.val[1], other.val[1]
        a2, b2 = self.val[2], other.val[2]

        v0 = a0 * b0
        v1 = a1 * b1
        v2 = a2 * b2

        return self.__class__([
            v0 + ((a1 + a2) * (b1 + b2) - (v1 + v2)) * self.non_residue,
            (((a0 + a1) * (b0 + b1)) - (v0 + v1)) + (v2 * self.non_residue),
            (((a0 + a2) * (b0 + b2)) - (v0 + v2)) + (v1),
        ])

    def __rmul__(self, other):
        return self * other

    def mul_by_fq2(self, fq2):
        return self.__class__([
            fq2 * self.val[0],
            fq2 * self.val[1],
            fq2 * self.val[2],
        ])

    def mul_scalar(self, base):
        res = Fq6.zero()
        rem = Fq(base)
        exp = self

        while rem.val != 0:
            if rem.val % 2 == 1:
                res = res + exp

            exp += exp
            rem.val= rem.val >> 1

        return res

    def inverse(self):
        a0 = self.val[0]
        a1 = self.val[1]
        a2 = self.val[2]

        t0 = a0 * a0
        t1 = a1 * a1
        t2 = a2 * a2
        t3 = a0 * a1
        t4 = a0 * a2
        t5 = a1 * a2

        c0 = t0 - (t5 * self.non_residue)
        c1 = (t2 * self.non_residue) - t3
        c2 = t1 - t4

        t6 = (a0 * c0) + (self.non_residue * ((a2 * c1) + (a1 * c2)))
        t6 = t6.inverse()

        return self.__class__([
            t6 * c0,
            t6 * c1,
            t6 * c2,
        ])

    def __div__(self, other):
        return self * other.inverse()

    def __truediv__(self, other):
        return self.__div__(other)

    def __repr__(self):
        return repr(self.val)

    def frobenius_map(self, power):
        assert(self.frobenius_coeffs_c1 is not None)
        assert(self.frobenius_coeffs_c2 is not None)
        return self.__class__([self.val[0].frobenius_map(power),
                               self.frobenius_coeffs_c1[power % 6] * self.val[1].frobenius_map(power),
                               self.frobenius_coeffs_c2[power % 6] * self.val[2].frobenius_map(power)])