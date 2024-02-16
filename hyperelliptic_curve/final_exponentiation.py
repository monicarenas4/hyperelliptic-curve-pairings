from _utils import frobenius_power, hamming_weight, NAf_hamming_weight, NAF
from math import log2, floor, ceil
from sage.all import Integer


def final_exponentiation_cp8(miller_fun, U, Fp, NAF_rep=False):
    """
    :param miller_fun: output of the miller loop
    :param U: 4-elements vector
    :param Fp:
    :param NAF_rep:
    :return: final pairing output
    """
    u, u0, lx, ly = U[0], U[1], U[2], U[3]

    miller_fun = miller_fun.frobenius(4) / miller_fun
    fp = miller_fun.frobenius(1)
    fu1 = miller_fun ** u
    f1 = fp * fu1
    fp2 = f1.frobenius(2)
    f1u1 = f1 ** u
    f1u2 = f1u1 ** u
    f2 = fp2 * f1u2
    f22 = f2 ** 2
    y1 = f22 ** u0
    y2 = y1 ** u0
    y3 = y2 ** u0
    y4 = y3 ** u0
    y32 = y3 ** 2
    y34 = y32 ** 2
    y38 = y34 ** 2
    y42 = y4 ** 2
    y44 = y42 ** 2
    y48 = y44 ** 2
    y416 = y48 ** 2
    y4162 = y416 ** 2
    y4164 = y4162 ** 2
    y4168 = y4164 ** 2
    y41616 = y4168 ** 2
    y41624 = y41616 * y4168
    y41625 = y41624 * y416
    y382 = y38 ** 2
    y384 = y382 ** 2
    y385 = y384 * y38
    x1 = y4168 * f2
    x2 = x1 ** lx
    x3 = x2 * y384
    x4 = x3 ** lx
    z1 = y41625 * y385
    z2 = z1 * f22
    N = miller_fun * x4 * z2
    y22 = y2 ** 2
    y23 = y22 * y2
    y22lx = y22 ** lx
    m = y22lx * y2 * y1
    m2 = m ** 2
    m4 = m2 ** 2
    M = m4 * y23
    t0 = N / M

    mk, sk, fk, sk_cyclo, ik = 27, 18, 6, 12, 69
    exp_u_n, exp_u0_n, exp_ulx_n, mult, sq, inv, frob_power = 3, 4, 3, 17, 17, 2, 3

    if not NAF_rep:
        exp_u = (floor(log2(u)) - 1) * sk_cyclo + (hamming_weight(Integer(u).digits(2)) - 1) * mk
        exp_u0 = (floor(log2(u0)) - 1) * sk_cyclo + (hamming_weight(Integer(u0).digits(2)) - 1) * mk
        exp_lx = (floor(log2(lx)) - 1) * sk_cyclo + (hamming_weight(Integer(lx).digits(2)) - 1) * mk
    else:
        exp_u = (ceil(log2(u)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(u)) - 1) * mk
        exp_u0 = (ceil(log2(u0)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(u0)) - 1) * mk
        exp_lx = (floor(log2(lx)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(lx)) - 1) * mk

    # count over Fp
    mult, sq, inv, frob = (mult * mk), (sq * sk), (inv * ik), (frob_power * fk)
    total_ops = (exp_u_n * exp_u) + (exp_u0_n * exp_u0) + (exp_ulx_n * exp_lx) + mult + sq + inv + frob

    return t0, exp_u, exp_u0, exp_lx, mult, sq, inv, frob, total_ops


def final_exponentiation_k16(f, U, W, k=16, NAF_rep: bool = False):
    """
    :param f:
    :param U:
    :param W:
    :param k:
    :param NAF_rep:
    :return:
    """
    u, um = U[0], U[1]
    f = frobenius_power(f, k, W, 8) / f
    f1 = f
    f2 = f1 ** 2
    f4 = f2 ** 2
    f8 = f4 ** 2
    fum1 = f ** um
    fu1 = fum1 * f1
    fum2 = fum1 ** um
    f2u1 = fu1 ** 2
    f4u1 = f2u1 ** 2
    fup2 = fum2 * f4u1
    g1 = fum2 ** u
    g2 = g1 ** u
    g3 = g2 ** u
    g4 = g3 ** u
    g5 = g4 ** 2
    fl0 = fup2 * g5
    fl0u = fl0 ** u
    fl0u2 = fl0u ** u
    fl0u3 = fl0u2 ** u
    fn = f8 * fl0u3
    fnu1 = fn ** u
    fnu2 = fnu1 ** u
    fnu3 = fnu2 ** u
    fnu4 = fnu3 ** u
    N0 = fl0
    N1 = frobenius_power(fn, k, W, 1)
    N2 = frobenius_power(fnu3, k, W, 2)
    N3 = frobenius_power(fl0u, k, W, 3)
    N4 = frobenius_power(fnu1, k, W, 4)
    N5 = frobenius_power(fnu4, k, W, 5)
    N6 = frobenius_power(fl0u2, k, W, 6)
    N7 = frobenius_power(fnu2, k, W, 7)

    N, M = (N0 * N2 * N3 * N5 * N6), (N1 * N4 * N7)
    t0 = N / M

    mk, sk, fk, sk_cyclo, ik = 81, 54, 14, 36, 159
    exp_u_n, exp_um_n, mult, sq, inv, frob_power = 11, 2, 12, 6, 2, 8

    if not NAF_rep:
        exp_u = (floor(log2(u)) - 1) * sk_cyclo + (hamming_weight(Integer(u).digits(2)) - 1) * mk
        exp_um = (floor(log2(um)) - 1) * sk_cyclo + (hamming_weight(Integer(um).digits(2)) - 1) * mk
    else:
        exp_u = (ceil(log2(u)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(u)) - 1) * mk
        exp_um = (ceil(log2(um)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(um)) - 1) * mk

    # count over Fp
    mult, sq, inv, frob = (mult * mk), (sq * sk), (inv * ik), (frob_power * fk)
    total_ops = (exp_u_n * exp_u) + (exp_um_n * exp_um) + mult + sq + inv + frob

    return t0, exp_u, exp_um, mult, sq, inv, frob, total_ops


def final_exponentiation_new_k16(f, U, W, k=16, NAF_rep: bool = False):
    """
    :param f:
    :param U:
    :param W:
    :param k:
    :param NAF_rep:
    :return:
    """
    u, up = U[0], U[1]
    t0 = 1
    f = frobenius_power(f, k, W, 8) / f
    f1 = f
    f2 = f1 ** 2
    f4 = f2 ** 2
    f8 = f4 ** 2
    fup1 = f ** up
    fup12 = fup1 ** up
    gu1 = fup12 ** u
    gu2 = gu1 ** u
    gu3 = gu2 ** u
    gu4 = gu3 ** u
    gu5 = gu4 ** u
    gu6 = gu5 ** u
    y1 = gu6
    y2 = frobenius_power(gu5, k, W, 8)  # 1 / gu5
    y3 = y1 * y2
    y4 = y3 ** 2
    y5 = y4 ** 2
    y6 = y5 * gu4
    y7 = y6 ** 2
    fl0 = f8 * fup12 * y7
    fl0u1 = fl0 ** u
    fl0u2 = fl0u1 ** u
    fl0u3 = fl0u2 ** u
    fl1 = f8 * fl0u3
    hu1 = fl1 ** u
    hu2 = hu1 ** u
    hu3 = hu2 ** u
    hu4 = hu3 ** u
    N0 = fl0
    N1 = frobenius_power(fl1, k, W, 1)
    N2 = frobenius_power(hu3, k, W, 2)
    N3 = frobenius_power(fl0u1, k, W, 3)
    N4 = frobenius_power(hu1, k, W, 4)
    N5 = frobenius_power(hu4, k, W, 5)
    N6 = frobenius_power(fl0u2, k, W, 6)
    N7 = frobenius_power(hu2, k, W, 7)
    N = N0 * N2 * N3 * N5 * N6
    M = N1 * N4 * N7
    t1 = frobenius_power(M, k, W, 8)
    t0 = N * t1  # / M

    mk, sk, fk, sk_cyclo, ik = 81, 54, 14, 36, 159
    exp_u_n, exp_up_n, mult, sq, inv, frob_power = 13, 2, 13, 6, 1, 8

    if not NAF_rep:
        exp_u = (floor(log2(u)) - 1) * sk_cyclo + (hamming_weight(Integer(u).digits(2)) - 1) * mk
        exp_up = (floor(log2(up)) - 1) * sk_cyclo + (hamming_weight(Integer(up).digits(2)) - 1) * mk
    else:
        exp_u = (ceil(log2(u)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(u)) - 1) * mk
        exp_up = (ceil(log2(up)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(up)) - 1) * mk

    # count over Fp
    mult, sq, inv, frob = (mult * mk), (sq * sk), (inv * ik), (frob_power * fk)
    total_ops = (exp_u_n * exp_u) + (exp_up_n * exp_up) + mult + sq + inv + frob

    return t0, exp_u, exp_up, mult, sq, inv, frob, total_ops


def final_exponentiation_k24(f, U, W, k=24):
    """
    :param f:
    :param U:
    :param W:
    :param k:
    :return:
    """
    p = 127986680668455487902583272776597527184931623516192302056052221652345318942213021208502589516925179965136000549331214706623410587603554498076547338902102067893684380712411361
    r = 63995560145816128259234865857573280107224142930766262992242874861039654232684076725049076391602255284139227431215841

    pow = (p ** k - 1) // r
    t0 = f ** pow

    exp_u, exp_up, mult, sq, inv, frob_power, total = 1, 1, 1, 1, 1, 1, 1

    return t0, exp_u, exp_up, mult, sq, inv, frob_power, total
