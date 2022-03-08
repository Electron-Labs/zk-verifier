from field.fq import Fq
from field.fq2 import Fq2

class AteG1PreComp:
    def __init__(self, x, y):
        assert(isinstance(x, Fq) and isinstance(y, Fq))
        self.px = x
        self.py = y

    def __eq__(self, other):
        if not isinstance(other, AteG1PreComp):
            return False

        return self.px == other.px and self.py == other.py

    def __repr__(self):
        return f"px: {self.px} py: {self.py}"

class AteEllCoeffs:
    def __init__(self, ell0, ellvw, ellvv):
        assert(isinstance(ell0, Fq2) and
               isinstance(ellvw, Fq2) and
               isinstance(ellvv, Fq2))
        self.ell0 = ell0
        self.ellvw = ellvw
        self.ellvv = ellvv

    def __eq__(self, other):
        if not isinstance(other == AteEllCoeffs):
            return False

        return self.ell0 == other.ell0 and self.ellvw == other.ellvw and self.ellvv == other.ellvv

    def __repr__(self):
        return f"ell0: {self.ell0}, ellvw: {self.ellvw}, ellvv: {self.ellvv}"

class AteG2PreComp:
    def __init__(self, x, y):
        assert(isinstance(x, Fq2) and isinstance(y, Fq2))
        self.qx = x
        self.qy = y
        self.coeffs = []

    def __eq__(self, other):
        if not isinstance(other, AteG2PreComp):
            return False

        return self.qx == other.qx and self.qy == other.qy and self.coeffs == other.coeffs

    def __repr__(self):
        return f"qx: {self.qx}, qy: {self.qy}, coeffs: {self.coeffs}"