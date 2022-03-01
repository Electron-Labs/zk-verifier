from field.fq import Fq
from field.fq2 import Fq2
from field.fq6 import Fq6

class G2:
    field_modulus = None
    non_residue = None

    def __init__(self, val):
        if not self.field_modulus:
            raise AttributeError("Field modulus is not set")

        if not self.non_residue:
            raise AttributeError("Quadratic non residue is not set")

        Fq.field_modulus = self.field_modulus
        Fq2.non_residue = Fq(self.non_residue)
        
        self.val = val
        if len(val) == 2:
            self.val.append(Fq2.one())

    
    @classmethod
    def zero(cls):
        return cls([Fq2.zero(), Fq2.one(), Fq2.zero()])

    def is_zero(self):
        return self.val[2] == Fq2.zero()

    def __add__(self, other):
        if self.is_zero():
            return other

        if other.is_zero():
            return self

        x1 = self.val[0]
        y1 = self.val[1]
        z1 = self.val[2]
        x2 = other.val[0]
        y2 = other.val[1]
        z2 = other.val[2]

        z1_2 = z1 * z1
        z2_2 = z2 * z2

        u1 = x1 * z2_2
        u2 = x2 * z1_2

        t0 = z2 * z2_2
        s1 = y1 * t0

        t1 = z1 * z1_2
        s2 = y2 * t1

        h = u2 - u1
        t2 = h + h
        i = t2 * t2
        j = h * i
        t3 = s2 - s1
        r = t3 + t3
        v = u1 * i
        t4 = r * r
        t5 = v + v
        t6 = t4 - j
        x3 = t6 - t5
        t7 = v - x3
        t8 = s1 * j
        t9 = t8 + t8
        t10 = r * t7

        y3 = t10 - t9

        t11 = z1 + z2
        t12 = t11 * t11
        t13 = t12 - z1_2
        t14 = t13 - z2_2
        z3 = t14 * h

        return self.__class__([x3, y3, z3])


    def __neg__(self):
        return self.__class__([self.val[0], -self.val[1], self.val[2]])
    
    def __sub__(self, other):
        return self + (-other)
    
    def double(self):
        if self.is_zero():
            return self

        x = self.val[0]
        y = self.val[1]
        z = self.val[2]

        a = x * x
        b = y * y
        c = b * b

        t0 = x + b
        t1 = t0 * t0
        t2 = t1 - a
        t3 = t2 - c

        d = t3 + t3
        e = a + a + a
        f = e * e

        t4 = d + d
        x3 = f - t4

        t5 = d - x3
        c_2 = c + c
        c_4 = c_2 + c_2
        t6 = c_4 + c_4
        t7 = e * t5
        y3 = t7 - t6

        t8 = y * z
        z3 = t8 + t8

        return self.__class__([x3, y3, z3])
    
    def __eq__(self, other):
        if self.is_zero():
            return other.is_zero()
    
        if other.is_zero():
            return self.is_zero()
    
        z1_2 = self.val[2] * self.val[2]
        z2_2 = other.val[2] * other.val[2]
    
        u1 = self.val[0] * z2_2
        u2 = other.val[0] * z1_2
    
        z1_3 = self.val[2] * z1_2
        z2_3 = other.val[2] * z2_2
    
        s1 = self.val[1] * z2_3
        s2 = other.val[1] * z1_3
    
        return u1 == u2 and s1 == s2
    
    def affine(self):
        if self.is_zero():
            return self.zero()
    
        zinv = Fq2.one() / self.val[2]
        zinv2 = zinv * zinv
        x = self.val[0] * zinv2
    
        zinv3 = zinv2 * zinv
        y = self.val[1] * zinv3
    
        return G2([x, y, Fq2.one()])
    
    def mul_scalar(self, base):
        q = G2([Fq2.zero(), Fq2.zero(), Fq2.zero()])
        d = Fq(base)
        r = self

        for i in range(d.bit_length(), -1, -1):
            q = q.double()
            if d.bit(i):
                q += r
    
        return q
