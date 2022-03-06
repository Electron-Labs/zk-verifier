from field.utils import mod_inverse

class Fq:

    field_modulus = None 

    def __init__(self, val):
        # Check if field_modulus is set or not
        if not self.field_modulus:
            raise AttributeError("Field Modulus hasn't been specified")

        if isinstance(val, Fq):
            self.val = val.val
        elif isinstance(val, int):
            self.val = val % self.field_modulus
        else:
            raise TypeError("Expected an int or FQ object, but got object of type {}".format(type(val)))

    @classmethod
    def zero(cls):
        return cls(0)

    @classmethod
    def one(cls):
        return cls(1)

    def __add__(self, other):
        on = other.val if isinstance(other, Fq) else Fq(other)
        return Fq((self.val + on) % self.field_modulus)

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        on = other.val if isinstance(other, Fq) else Fq(other)
        return Fq((self.val * on) % self.field_modulus)

    def __rmul__(self, other):
        return self * other

    def __sub__(self, other):
        on = other.val if isinstance(other, Fq) else Fq(other)
        return  Fq((self.val - on) % self.field_modulus)

    def __rsub__(self, other):
        on = other.val if isinstance(other, Fq) else Fq(other)
        return  Fq((on - self.val) % self.field_modulus)

    def __div__(self, other):
        on = other.val if isinstance(other, Fq) else Fq(other)
        return Fq((self.val * mod_inverse(on, self.field_modulus)) % self.field_modulus)

    def __truediv__(self, other):
        return self.__div__(other)

    def __rdiv__(self, other):
        on = other.val if isinstance(other, Fq) else Fq(other)
        return Fq((mod_inverse(self.val, self.field_modulus) * on) % self.field_modulus)

    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    def __pow__(self, other):
        if other == 0:
            return Fq.one()
        elif other == 1:
            return Fq(self.val)
        elif other % 2 == 0:
            return (self * self) ** (other // 2)
        else:
            return ((self * self) ** int(other // 2)) * self

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return Fq(-self.val)

    def __repr__(self):
        return repr(self.val)

    def __eq__(self, other):
        if isinstance(other, Fq):
            return self.val == other.val
        else:
            return self.val == other

    def bit_length(self):
        return self.val.bit_length()

    def bit(self, i):
        if self.val & (1 << i):
            return True
        else:
            return False