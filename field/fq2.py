from field.fq import Fq

class Fq2:
    non_residue = None
    degree = 2

    def __init__(self, val):
        if not self.non_residue:
            raise AttributeError("Quadratic non-residue hasn't been specified")

        assert(isinstance(self.non_residue, Fq))
        assert(isinstance(val, list))

        if isinstance(val, Fq2):
            self.val = val
        else:
            self.val = [Fq(v) for v in val]

        assert(len(self.val) == self.degree)

    @classmethod
    def zero(cls):
        return cls([Fq.zero(), Fq.zero()])

    @classmethod
    def one(cls):
        return cls([Fq.one(), Fq.zero()])

    def __eq__(self, other):
        assert isinstance(other, self.__class__)
        for c1, c2 in zip(self.val, other.val):
            if c1 != c2:
                return False
        return True

    def __add__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x + y for x, y in zip(self.val, other.val)])

    def __sub__(self, other):
        assert isinstance(other, self.__class__)
        return self.__class__([x - y for x, y in zip(self.val, other.val)])

    def __mul__(self, other):
        # Multiplication and Squaring on Pairing-Friendly.pdf; Section 3 (Karatsuba)
        # https://pdfs.semanticscholar.org/3e01/de88d7428076b2547b60072088507d881bf1.pdf

        a0 = self.val[0]
        a1 = self.val[1]
        b0 = other.val[0]
        b1 = other.val[1]

        v0 = self.val[0] * other.val[0]
        v1 = self.val[1] * other.val[1]

        val1 = v0 + v1 * self.non_residue
        val2 = (a0 + a1) * (b0 + b1) - (v0 + v1)

        return self.__class__([val1, val2])

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return type(self)([-c for c in self.val])

    def square(self):
        a0 = self.val[0]
        a1 = self.val[1]

        ab = a0 * a1

        return self.__class__([
            ((a0 + a1) * (a0 + self.non_residue * a1)) - (ab + self.non_residue * ab),
            ab + ab,
        ])

    def inverse(self):
        # High-Speed Software Implementation of the Optimal Ate Pairing over Barreto–Naehrig Curves .pdf
        # https://eprint.iacr.org/2010/354.pdf , algorithm 8
        t0 = self.val[0] * self.val[0]
        t1 = self.val[1] * self.val[1]
        t2 = t0 - (t1 * self.non_residue)
        t3 = Fq.one() / t2
        return self.__class__([(self.val[0] * t3), -(self.val[1] * t3)])

    def __div__(self, other):
        return self * other.inverse()

    def __repr__(self):
        return repr(self.val)

    def mul_scalar(self, base):
        q = Fq2.zero()
        d = Fq(base)
        r = self

        found_one = False

        for i in range(d.bit_length(), -1, -1):
            if found_one:
                q += q
            if d.bit(i):
                found_one = True
                q += r

        return q
