from xmlrpc.server import resolve_dotted_attribute
import pytest
from field.fq import Fq
from field.fq2 import Fq2
from field.fq6 import Fq6
from field.fq12 import Fq12


def test_fq12_add():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    a1 = Fq6([Fq2([10, 2]), Fq2([1, 23]), Fq2([12, 24])])
    b1 = Fq6([Fq2([2, 21]), Fq2([54, 21]), Fq2([110, 52])])
    c = Fq12([a, b])
    c1 = Fq12([a1, b1])
    res_a = Fq6([Fq2([11, 4]), Fq2([2, 25]), Fq2([13, 26])])
    res_b = Fq6([Fq2([4, 23]), Fq2([59, 23]), Fq2([120, 54])])
    res = Fq12([res_a, res_b])
    d = c + c1
    assert(d == res)


def test_fq12_sub():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    a1 = Fq6([Fq2([10, 2]), Fq2([1, 23]), Fq2([12, 24])])
    b1 = Fq6([Fq2([2, 21]), Fq2([54, 21]), Fq2([110, 52])])
    c = Fq12([a, b])
    c1 = Fq12([a1, b1])
    res_a = Fq6([Fq2([9, 0]), Fq2([0, 21]), Fq2([11, 22])])
    res_b = Fq6([Fq2([0, 19]), Fq2([49, 19]), Fq2([100, 50])])
    res = Fq12([res_a, res_b])
    d = c1 - c
    assert(d == res)


def test_fq12_mul():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    a1 = Fq6([Fq2([10, 2]), Fq2([1, 23]), Fq2([12, 24])])
    b1 = Fq6([Fq2([2, 21]), Fq2([54, 21]), Fq2([110, 52])])
    c = Fq12([a, b])
    c1 = Fq12([a1, b1])
    res_a = Fq6([Fq2([1351, 7679]), Fq2([7249, 8615]), Fq2([8183, 8010])])
    res_b = Fq6([Fq2([21888242871839275222246405745257275088696311157297823662689037894645226207728, 7036]), Fq2([140, 5134]), Fq2([9, 655])])
    res = Fq12([res_a, res_b])
    d = c * c1
    assert(d == res)


def test_fq12_inverse():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    c = Fq12([a, b])
    c_inv = c.inverse()
    one = Fq12.one()
    assert(c * c_inv == one)


def test_fq12_div():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    a1 = Fq6([Fq2([10, 2]), Fq2([1, 23]), Fq2([12, 24])])
    b1 = Fq6([Fq2([2, 21]), Fq2([54, 21]), Fq2([110, 52])])
    c = Fq12([a, b])
    c1 = Fq12([a1, b1])
    d = c * c1
    e = d / c1
    assert(e == c)


def test_fq12_pow():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    c = Fq12([a, b])
    c_2 = c ** 2
    res = c * c
    assert(c_2 == res)
    c_3 = c ** 3
    res *= c
    assert(c_3 == res)


def test_fq12_mul_scalar():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    c = Fq12([a, b])
    d = c.mul_scalar(4)
    res_a = Fq6([Fq2([4, 8]), Fq2([4, 8]), Fq2([4, 8])])
    res_b = Fq6([Fq2([8, 8]), Fq2([20, 8]), Fq2([40, 8])])
    res = Fq12([res_a, res_b])
    assert(d == res)


def test_fq12_frobenius_map():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    Fq2.frobenius_coeffs_c1 = [Fq(1),
                               Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)]
    Fq6.frobenius_coeffs_c1 = [Fq2([1, 0]),
                               Fq2([21575463638280843010398324269430826099269044274347216827212613867836435027261, 10307601595873709700152284273816112264069230130616436755625194854815875713954]),
                               Fq2([21888242871839275220042445260109153167277707414472061641714758635765020556616, 0]),
                               Fq2([3772000881919853776433695186713858239009073593817195771773381919316419345261, 2236595495967245188281701248203181795121068902605861227855261137820944008926]),
                               Fq2([2203960485148121921418603742825762020974279258880205651966, 0]),
                               Fq2([18429021223477853657660792034369865839114504446431234726392080002137598044644, 9344045779998320333812420223237981029506012124075525679208581902008406485703]),
                              ]
    Fq6.frobenius_coeffs_c2 = [Fq2([1, 0]),
                               Fq2([2581911344467009335267311115468803099551665605076196740867805258568234346338, 19937756971775647987995932169929341994314640652964949448313374472400716661030]),
                               Fq2([2203960485148121921418603742825762020974279258880205651966, 0]),
                               Fq2([5324479202449903542726783395506214481928257762400643279780343368557297135718, 16208900380737693084919495127334387981393726419856888799917914180988844123039]),
                               Fq2([21888242871839275220042445260109153167277707414472061641714758635765020556616, 0]),
                               Fq2([13981852324922362344252311234282257507216387789820983642040889267519694726527, 7629828391165209371577384193250820201684255241773809077146787135900891633097]),
                              ]
    Fq12.frobenius_coeffs_c1 = [Fq2([1, 0]),
                                Fq2([8376118865763821496583973867626364092589906065868298776909617916018768340080, 16469823323077808223889137241176536799009286646108169935659301613961712198316]),
                                Fq2([21888242871839275220042445260109153167277707414472061641714758635765020556617, 0]),
                                Fq2([11697423496358154304825782922584725312912383441159505038794027105778954184319, 303847389135065887422783454877609941456349188919719272345083954437860409601]),
                                Fq2([21888242871839275220042445260109153167277707414472061641714758635765020556616, 0]),
                                Fq2([3321304630594332808241809054958361220322477375291206261884409189760185844239, 5722266937896532885780051958958348231143373700109372999374820235121374419868]),
                                Fq2([21888242871839275222246405745257275088696311157297823662689037894645226208582, 0]),
                                Fq2([13512124006075453725662431877630910996106405091429524885779419978626457868503, 5418419548761466998357268504080738289687024511189653727029736280683514010267]),
                                Fq2([2203960485148121921418603742825762020974279258880205651966, 0]),
                                Fq2([10190819375481120917420622822672549775783927716138318623895010788866272024264, 21584395482704209334823622290379665147239961968378104390343953940207365798982]),
                                Fq2([2203960485148121921418603742825762020974279258880205651967, 0]),
                                Fq2([18566938241244942414004596690298913868373833782006617400804628704885040364344, 16165975933942742336466353786298926857552937457188450663314217659523851788715]),
                               ]
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    c = Fq12([a, b])
    c = c.frobenius_map(3)

    res_a = Fq6([Fq2([1, 21888242871839275222246405745257275088696311157297823662689037894645226208581]),
                 Fq2([8245191873854344152997097683120221829251211399028918227483904194958307363113, 16580836603966812857660716620032740405799232872269293346997535193833331526987]),
                 Fq2([15854037092086014490319367904917715356019399444816597216927133835889759173213, 5559941975837885999465928336321959017537210895055602240357227443874249851603])])
    res_b = Fq6([Fq2([2114298899147165162250727009667395420041154102860624959589184225788402979257, 20989333529232373609686812555100319434480553810116075792480189486608264867730]),
                 Fq2([768729425043478242131547316821264100817738211278702545271283937646543715007, 10002154396670567017907573318723036179714327978289766191776806537852539809676]),
                 Fq2([647599598232173817913279931440377159034147622558037926956600255118455891821, 19834765063218638114715815552011489526788980899606756647125693726673574182815])])
    res = Fq12([res_a, res_b])
    assert(c == res)


def test_fq12_cyclomatic_squared():
    Fq.field_modulus = 21888242871839275222246405745257275088696311157297823662689037894645226208583
    Fq2.non_residue = Fq(21888242871839275222246405745257275088696311157297823662689037894645226208582)
    Fq6.non_residue = Fq2([9, 1])
    Fq12.non_residue = Fq2([9, 1])
    a = Fq6([Fq2([1, 2]), Fq2([1, 2]), Fq2([1, 2])])
    b = Fq6([Fq2([2, 2]), Fq2([5, 2]), Fq2([10, 2])])
    c = Fq12([a, b])
    d = c.cyclotomic_square()
    res_a = Fq6([Fq2([496, 611]),
                 Fq2([21888242871839275222246405745257275088696311157297823662689037894645226208488, 119]),
                 Fq2([2461, 1376])])
    res_b = Fq6([Fq2([196, 1228]),
                 Fq2([16, 76]),
                 Fq2([8, 40])])
    res = Fq12([res_a, res_b])

    assert(d == res)