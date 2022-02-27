from field.utils import mod_inverse

class Fq:

    field_modulus = None 

    def __init__(self, val):
        # Check if field_modulus is set or not
        if not self.field_modulus:
            raise Exception("field_modulus is not set")

        self.val = val % self.field_modulus

    @classmethod
    def zero(cls):
        return cls(0)

    @classmethod
    def one(cls):
        return cls(1)

    def add(self, other):
        assert(isinstance(other, Fq))
        self.val = (self.val + other.val) % self.field_modulus

    def double(self):
        self.val = (self.val + self.val) % self.field_modulus

    def sub(self, other):
        assert(isinstance(other, Fq))
        self.val = (self.val - other.val) % self.field_modulus

    def mul(self, other):
        assert(isinstance(other, Fq))
        self.val = (self.val * other.val) % self.field_modulus

    def mul_scalar(self, e):
        self.mul(Fq(e))

    def neg(self):
        self.mul_scalar(-1)

    def inverse(self):
        self.val = mod_inverse(self.val, self.field_modulus)

    def div(self, other):
        assert(isinstance(other, Fq))
        copy = other
        copy.inverse()
        assert(isinstance(copy, Fq))
        self.mul(copy)

    def exp(self, exponent):
        pow(self.val, exponent, self.field_modulus)

    def square(self):
        self.mul(self)

    def eq(self, other):
        assert(isinstance(other, Fq))
        return (self.val == other.val)
