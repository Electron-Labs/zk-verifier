from bn128.ate import AteEllCoeffs, AteG1PreComp, AteG2PreComp
from field.fq import Fq
from field.fq2 import Fq2


def test_ate_g1_precomp_repr():
    Fq.field_modulus = 17
    ate_g1 = AteG1PreComp(Fq(1), Fq(2))
    assert(repr(ate_g1) == "px: 1 py: 2")


def test_ate_ell_coeffs_repr():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    ell_coeffs = AteEllCoeffs(Fq2([1, 2]), Fq2([4, 5]), Fq2([9, 7]))
    assert(repr(ell_coeffs) == "ell0: [1, 2], ellvw: [4, 5], ellvv: [9, 7]")


def test_ate_g2_precomp_repr():
    Fq.field_modulus = 17
    Fq2.non_residue = Fq(16)
    ell_coeffs = AteEllCoeffs(Fq2([1, 2]), Fq2([4, 5]), Fq2([9, 7]))
    ate_g2 = AteG2PreComp(Fq2([11, 22]), Fq2([41, 35]))
    ate_g2.coeffs.append(ell_coeffs)
    ate_g2.coeffs.append(ell_coeffs)
    assert(repr(ate_g2) == "qx: [11, 5], qy: [7, 1], coeffs: [ell0: [1, 2], ellvw: [4, 5], ellvv: [9, 7], ell0: [1, 2], ellvw: [4, 5], ellvv: [9, 7]]")