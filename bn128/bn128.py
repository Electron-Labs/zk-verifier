from field.fq import Fq
from field.fq2 import Fq2
from field.fq6 import Fq6
from field.fq12 import Fq12
from bn128.g1 import G1
from bn128.g2 import G2
from bn128.ate import *


class BN128:
    def __init__(self):
        # Parameters for base field (Fq)
        Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583

        # Parameters for twist field (Fq2)
        Fq2.non_residue = Fq(
            21888242871839275222246405745257275088696311157297823662689037894645226208582
        )
        Fq2.frobenius_coeffs_c1 = [
            Fq(1),
            Fq(
                21888242871839275222246405745257275088696311157297823662689037894645226208582
            ),
        ]

        # Parameters for Fq6
        Fq6.non_residue = Fq2([9, 1])
        Fq6.frobenius_coeffs_c1 = [
            Fq2([1, 0]),
            Fq2(
                [
                    21575463638280843010398324269430826099269044274347216827212613867836435027261,
                    10307601595873709700152284273816112264069230130616436755625194854815875713954,
                ]
            ),
            Fq2(
                [
                    21888242871839275220042445260109153167277707414472061641714758635765020556616,
                    0,
                ]
            ),
            Fq2(
                [
                    3772000881919853776433695186713858239009073593817195771773381919316419345261,
                    2236595495967245188281701248203181795121068902605861227855261137820944008926,
                ]
            ),
            Fq2([2203960485148121921418603742825762020974279258880205651966, 0]),
            Fq2(
                [
                    18429021223477853657660792034369865839114504446431234726392080002137598044644,
                    9344045779998320333812420223237981029506012124075525679208581902008406485703,
                ]
            ),
        ]
        Fq6.frobenius_coeffs_c2 = [
            Fq2([1, 0]),
            Fq2(
                [
                    2581911344467009335267311115468803099551665605076196740867805258568234346338,
                    19937756971775647987995932169929341994314640652964949448313374472400716661030,
                ]
            ),
            Fq2([2203960485148121921418603742825762020974279258880205651966, 0]),
            Fq2(
                [
                    5324479202449903542726783395506214481928257762400643279780343368557297135718,
                    16208900380737693084919495127334387981393726419856888799917914180988844123039,
                ]
            ),
            Fq2(
                [
                    21888242871839275220042445260109153167277707414472061641714758635765020556616,
                    0,
                ]
            ),
            Fq2(
                [
                    13981852324922362344252311234282257507216387789820983642040889267519694726527,
                    7629828391165209371577384193250820201684255241773809077146787135900891633097,
                ]
            ),
        ]

        # Parameters for Fq12
        Fq12.non_residue = Fq2([9, 1])
        Fq12.frobenius_coeffs_c1 = [
            Fq2([1, 0]),
            Fq2(
                [
                    8376118865763821496583973867626364092589906065868298776909617916018768340080,
                    16469823323077808223889137241176536799009286646108169935659301613961712198316,
                ]
            ),
            Fq2(
                [
                    21888242871839275220042445260109153167277707414472061641714758635765020556617,
                    0,
                ]
            ),
            Fq2(
                [
                    11697423496358154304825782922584725312912383441159505038794027105778954184319,
                    303847389135065887422783454877609941456349188919719272345083954437860409601,
                ]
            ),
            Fq2(
                [
                    21888242871839275220042445260109153167277707414472061641714758635765020556616,
                    0,
                ]
            ),
            Fq2(
                [
                    3321304630594332808241809054958361220322477375291206261884409189760185844239,
                    5722266937896532885780051958958348231143373700109372999374820235121374419868,
                ]
            ),
            Fq2(
                [
                    21888242871839275222246405745257275088696311157297823662689037894645226208582,
                    0,
                ]
            ),
            Fq2(
                [
                    13512124006075453725662431877630910996106405091429524885779419978626457868503,
                    5418419548761466998357268504080738289687024511189653727029736280683514010267,
                ]
            ),
            Fq2([2203960485148121921418603742825762020974279258880205651966, 0]),
            Fq2(
                [
                    10190819375481120917420622822672549775783927716138318623895010788866272024264,
                    21584395482704209334823622290379665147239961968378104390343953940207365798982,
                ]
            ),
            Fq2([2203960485148121921418603742825762020974279258880205651967, 0]),
            Fq2(
                [
                    18566938241244942414004596690298913868373833782006617400804628704885040364344,
                    16165975933942742336466353786298926857552937457188450663314217659523851788715,
                ]
            ),
        ]

        G1.field_modulus = Fq.field_modulus
        G2.field_modulus = Fq.field_modulus
        G2.non_residue = Fq2.non_residue

        # y ** 2 = x ** 3 + b (Ellipitic Curve Equation)
        self.b = Fq(3)

        # Twist curve parameters
        self.twist = Fq2([9, 1])
        self.twist_b = self.twist.inverse().mul_scalar(self.b)
        self.two_inv = Fq.one() / Fq(2)
        self.twist_mul_by_x = Fq2(
            [
                21575463638280843010398324269430826099269044274347216827212613867836435027261,
                10307601595873709700152284273816112264069230130616436755625194854815875713954,
            ]
        )
        self.twist_mul_by_y = Fq2(
            [
                2821565182194536844548159561693502659359617185244120367078079554186484126554,
                3505843767911556378687030309984248845540243509899259641013678093033130930403,
            ]
        )

        # Ate pairing constans
        self.final_exp = 552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480
        self.final_exp_z = 4965661367192848881
        self.ate_loop_count = Fq(29793968203157093288)
        self.ate_is_loop_count_neg = False
        self.final_exp_is_z_neg = False

        # Generators
        self.g1_zero = G1.zero()
        self.g1 = G1([1, 2])

        self.g2_zero = G2.zero()
        self.g2 = G2(
            [
                Fq2(
                    [
                        10857046999023057135944570762232829481370756359578518086990519993285655852781,
                        11559732032986387107991004021392285783925812861821192530917403151452391805634,
                    ]
                ),
                Fq2(
                    [
                        8495653923123431417604973247489272438418190587263600148770280649306958101930,
                        4082367875863433681332203403145435568316851327593401208105741076214120093531,
                    ]
                ),
            ]
        )

    def ate_precompute_g1(self, p):
        print("Entering precompute G1")
        assert(isinstance(p, G1))
        pcopy = p.affine()
        ate_g1_precomp = AteG1PreComp(pcopy.val[0], pcopy.val[1])
        print("End with precompute G1")
        return ate_g1_precomp

    def doubling_step_for_flipped_miller_loop(self, curr_g2):
        x = curr_g2.val[0]
        y = curr_g2.val[1]
        z = curr_g2.val[2]

        a = x * y
        a = a.mul_scalar(self.two_inv)
        b = y * y
        c = z * z
        d = c + c + c
        e = self.twist_b * d
        f = e + e + e
        g = (b + f)
        g = g.mul_scalar(self.two_inv)
        h = ((y + z) * (y + z)) - (b + c)
        i = e - b
        j = x * x
        e_sq = e * e

        res_x = a * (b - f)
        res_y = (g * g) - (e_sq + e_sq + e_sq)
        res_z = b * h

        ell0 = i * self.twist
        ellvw = -h
        ellvv = j + j + j

        ate_ell_coeff = AteEllCoeffs(ell0, ellvw, ellvv)

        return G2([res_x, res_y, res_z]), ate_ell_coeff

    def mul_by_q(self, q):
        assert(isinstance(q, G2))

        fmx = self.twist_mul_by_x * q.val[0].frobenius_map(1)
        fmy = self.twist_mul_by_y * q.val[1].frobenius_map(1)
        fmz = q.val[2].frobenius_map(1)

        return G2([fmx, fmy, fmz])

    def mixed_addition_step_for_flipped_miller(self, base_g2, curr_g2):
        x1 = curr_g2.val[0]
        y1 = curr_g2.val[1]
        z1 = curr_g2.val[2]
        x2 = base_g2.val[0]
        y2 = base_g2.val[1]

        d = x1 - (x2 * z1)
        e = y1 - (y2 * z1)
        f = d * d
        g = e * e
        h = d * f
        i = x1 * f
        j = (h + (z1 * g)) - (i + i)

        res_x = d * j
        res_y = (e * (i - j)) - (h * y1)
        res_z = z1 * h

        ell0 = self.twist * ((e * x2) - (d * y2))
        ellvv = -e
        ellvw = d

        ate_ell_coeff = AteEllCoeffs(ell0, ellvw, ellvv)

        return G2([res_x, res_y, res_z]), ate_ell_coeff

    def ate_precompute_g2(self, q):
        print("Entering precompute G2")
        assert(isinstance(q, G2))
        qcopy = q.affine()
        ate_g2_precomp = AteG2PreComp(qcopy.val[0], qcopy.val[1])

        r = G2([qcopy.val[0], qcopy.val[1], Fq2.one()])

        found_one = False

        for i in range(self.ate_loop_count.bit_length(), -1, -1):
            bit = self.ate_loop_count.bit(i)

            if not found_one:
                found_one |= bit
                continue

            r, c = self.doubling_step_for_flipped_miller_loop(r)
            ate_g2_precomp.coeffs.append(c)

            if bit:
                r, c = self.mixed_addition_step_for_flipped_miller(qcopy, r)
                found_one = True
                ate_g2_precomp.coeffs.append(c)

        q1 = self.mul_by_q(qcopy).affine()
        assert(q1.val[2] == Fq2.one())

        q2 = self.mul_by_q(q1).affine()
        assert(q2.val[2] == Fq2.one())

        if self.ate_is_loop_count_neg:
            r.val[1] = -r.val[1]

        q2.val[1] = -q2.val[1]

        r, c = self.mixed_addition_step_for_flipped_miller(q1, r)
        ate_g2_precomp.coeffs.append(c)

        r, c = self.mixed_addition_step_for_flipped_miller(q2, r)
        ate_g2_precomp.coeffs.append(c)

        print("Done with precompute G2")

        return ate_g2_precomp

    def ate_miller_loop(self, p, q):
        assert(isinstance(p, AteG1PreComp) and isinstance(q, AteG2PreComp))
        f = Fq12.one()
        found_one = False
        idx = 0
        c = q.coeffs[idx]

        for i in range(self.ate_loop_count.bit_length(), -1, -1):
            bit = self.ate_loop_count.bit(i)
            if not found_one:
                found_one |= bit
                continue
            c = q.coeffs[idx]
            assert(isinstance(c, AteEllCoeffs))
            idx += 1
            f = f.square()
            f = f.mul_by_024(c.ell0, c.ellvw.mul_scalar(p.py), c.ellvv.mul_scalar(p.px))

            if bit:
                c = q.coeffs[idx]
                f = f.mul_by_024(c.ell0, c.ellvw.mul_scalar(p.py), c.ellvv.mul_scalar(p.px))
                idx += 1

        if self.ate_is_loop_count_neg:
            f = f.inverse()

        c = q.coeffs[idx]
        idx += 1
        f = f.mul_by_024(c.ell0, c.ellvw.mul_scalar(p.py), c.ellvv.mul_scalar(p.px))

        c = q.coeffs[idx]
        idx += 1
        f = f.mul_by_024(c.ell0, c.ellvw.mul_scalar(p.py), c.ellvv.mul_scalar(p.px))

        return f

    def final_exp_first_chunk(self, elt):
        assert(isinstance(elt, Fq12))
        a = elt.unitary_inverse()
        b = elt.inverse()
        c = a * b
        d = c.frobenius_map(2)
        result = d * c

        return result

    def final_exp_by_neg_z(self, elt):
        assert(isinstance(elt, Fq12))
        result = elt.cyclotomic_exp(self.final_exp_z)
        if not self.final_exp_is_z_neg:
            result = result.unitary_inverse()

        return result

    def final_exp_last_chunk(self, elt):
        assert(isinstance(elt, Fq12))
        a = self.final_exp_by_neg_z(elt)
        b = a.cyclotomic_square()
        c = b.cyclotomic_square()
        d = c * b
        e = self.final_exp_by_neg_z(d)
        f = e.cyclotomic_square()
        g = self.final_exp_by_neg_z(f)
        h = d.unitary_inverse()
        i = g.unitary_inverse()
        j = i * e
        k = j * h
        l = k * b
        m = k * e
        n = m * elt
        o = l.frobenius_map(1)
        p = o * n
        q = k.frobenius_map(2)
        r = q * p
        s = elt.unitary_inverse()
        t = s * l
        u = t.frobenius_map(3)
        v = u * r

        result = v

        return result

    def final_exponentation(self, elt):
        assert(isinstance(elt, Fq12))
        a = self.final_exp_first_chunk(elt)
        result = self.final_exp_last_chunk(a)

        return result

    def pairing_check(self, p, q):
        assert(isinstance(p, G1) and isinstance(q, G2))
        print("Starting pairing check")
        precomp_g1 = self.ate_precompute_g1(p)
        assert(isinstance(precomp_g1, AteG1PreComp))
        precomp_g2 = self.ate_precompute_g2(q)
        assert(isinstance(precomp_g2, AteG2PreComp))
        result = self.ate_miller_loop(precomp_g1, precomp_g2)
        result = self.final_exponentation(result)
        print("Done with pairing check")
        return result