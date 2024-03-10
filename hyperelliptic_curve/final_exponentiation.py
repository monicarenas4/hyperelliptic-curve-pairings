from _utils import frobenius_power, hamming_weight, NAf_hamming_weight, NAF
from math import log2, floor, ceil
from sage.all import Integer


def final_exponentiation_cp8(miller_fun, U, Fp, NAF_rep=False):
    """
    :param miller_fun: Miller function, output of the miller loop
    :param U: 4-elements vector U = [u, u/4, lx, ly]
    :param Fp: base field
    :param NAF_rep: determine whether NAF or binary representation is used
    :return: final pairing output

    easy part: f^(p^4 - 1)
    cost: 1 mult. + 1 inv. + 1 Frobenius

    hard part: f^((p^4 + 1)/r)
    hard part formula:
    e = (p^4 + 1)/r = (lx^2 + 2*ly^2 + ly + 1/8)*u^4 + (lx + ly + 1/4)*u^3 - (lx + ly - 1/8)*u^2 - (ly + 1)*u + (lx^2 + 2*ly^2)
    cost: 3 exp_u + 4 exp_u0 + 3 exp_lx + 16 mult. + 17 sq. + 1 inv. + 2 Frobenius
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


def final_exponentiation_KT16(f, U, W, k=16, NAF_rep: bool = False):
    """
    :param f: Miller function, output of the miller loop
    :param U: 2-elements vector U = [u, u - 1]
    :param W: powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    :param k: embedding degree
    :param NAF_rep: determine whether NAF or binary representation is used
    :return: final pairing output

    easy part: f^(p^8 - 1)
    cost: 1 mult. + 1 inv. + 1 Frobenius

    hard part: f^((p^8 + 1)/r)
    hard part formula:
    e = (p^8 + 1)/r = l0 + l1*p + l2*p^2 + l3*p^3 + l4*p^4 + l5*p^5 + l6*p^6 + l7*p^7

    l0 = (u + 1)^2 + 2*u^4*(u - 1)^2
    l1 = - 8 - u^3*l0
    l2 = - u^3*l1
    l3 = u*l0
    l4 = u*l1
    l5 = u*l2
    l6 = u*l3
    l7 = u*l4
    cost: 11 exp_u + 2 exp_(u-1) + 11 mult. + 6 sq. + 1 inv. + 7 Frobenius
    """
    d = 8 # degre of twist
    u, um = U[0], U[1]
    f = frobenius_power(f, k, W, d, 8) / f
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
    N1 = frobenius_power(fn, k, W, d, 1)
    N2 = frobenius_power(fnu3, k, W, d, 2)
    N3 = frobenius_power(fl0u, k, W, d, 3)
    N4 = frobenius_power(fnu1, k, W, d, 4)
    N5 = frobenius_power(fnu4, k, W, d, 5)
    N6 = frobenius_power(fl0u2, k, W, d, 6)
    N7 = frobenius_power(fnu2, k, W, d, 7)

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


def final_exponentiation_New16(f, U, W, k=16, NAF_rep: bool = False):
    """
    :param f: Miller function, output of the miller loop
    :param U: 2-elements vector U = [u, u + 1]
    :param W: powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    :param k: embedding degree
    :param NAF_rep: determine whether NAF or binary representation is used
    :return: final pairing output

    easy part: f^(p^8 - 1)
    cost: 1 mult. + 1 inv. + 1 Frobenius

    hard part: f^((p^8 + 1)/r)
    hard part formula:
    e = (p^8 + 1)/r = l0 + l1*p + l2*p^2 + l3*p^3 + l4*p^4 + l5*p^5 + l6*p^6 + l7*p^7

    l0 = 8 + (u + 1)^2 * (1 + 2*u^4 - 8*u^5 + 8*u^6)
    l1 = - (8 + u^3*l0)
    l2 = u^3*l1
    l3 = u*l0
    l4 = u*l1
    l5 = -u^4*l1
    l6 = u^2*l2
    l7 = u^2*l1
    cost: 13 exp_u + 2 exp_(u+1) + 12 mult. + 6 sq. + 1 inv. + 7 Frobenius
    """
    d = 8 # degree of twist
    u, up = U[0], U[1]
    t0 = 1
    f = frobenius_power(f, k, W, d, 8) / f
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
    y2 = frobenius_power(gu5, k, W, d, 8)  # 1 / gu5
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
    N1 = frobenius_power(fl1, k, W, d, 1)
    N2 = frobenius_power(hu3, k, W, d, 2)
    N3 = frobenius_power(fl0u1, k, W, d, 3)
    N4 = frobenius_power(hu1, k, W, d, 4)
    N5 = frobenius_power(hu4, k, W, d, 5)
    N6 = frobenius_power(fl0u2, k, W, d, 6)
    N7 = frobenius_power(hu2, k, W, d, 7)
    N = N0 * N2 * N3 * N5 * N6
    M = N1 * N4 * N7
    t1 = frobenius_power(M, k, W, d, 8)
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


def final_exponentiation_New24(f, U, W, k=24, NAF_rep: bool = False):
    """
    :param f: Miller function, output of the miller loop
    :param U: 1-element vector U = [u]
    :param W: powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    :param k: embedding degree
    :return: final pairing output

    easy part: f^(p^12 - 1)*(p^4 + 1)
    cost: 2 mult. + 1 inv. + 2 Frobenius

    hard part: f^((p^8 - p^4 + 1)/r)
    hard part formula:
    e = (p^8 - p^4 + 1)/r = l0 + l1*p + l2*p^2 + l3*p^3 + l4*p^4 + l5*p^5 + l6*p^6 + l7*p^7

    l0 = u*l7
    l1 = -u*l0
    l2 = -u*l1
    l3 = -u*l2
    l4 = -1 + 2*u - 3*u^2 - 4*u^3 - 2*u^4
    l5 = 8 - u*l4
    l6 = -u*l5
    l7 = -u*l6
    cost: 11 exp_u + 13 mult. + 8 sq. + 1 inv. + 11 Frobenius
    """
    d = 8 # degree of twist

    u = U[0]

    # easy part: f^(p^12 - 1)*(p^4 + 1)
    fp12 = frobenius_power(f, k, W, d, 12)
    invf = 1 / f
    f = fp12 * invf
    f = f * frobenius_power(f, k, W, d, 4)

    # hard part: f^(p^8 - p^4 + 1)
    fu1 = f ** u
    fu2 = fu1 ** u
    fu3 = fu2 ** u
    fu4 = fu3 ** u
    t1 = fu1 ** 2
    t2 = fu2 ** 2
    t2 = t2 * fu2
    t3 = fu3 ** 2
    t3 = t3 ** 2
    t4 = fu4 ** 2
    g1 = frobenius_power(t1, k, W, d, 12)
    fl4 = f * g1 * t2 * t3 * t4

    f2 = f ** 2
    f4 = f2 ** 2
    f8 = f4 ** 2
    fl5 = fl4 ** u
    fl5 = fl5 * f8

    fl6 = fl5 ** u
    fl7 = fl6 ** u
    fl0 = fl7 ** u
    fl1 = fl0 ** u
    fl2 = fl1 ** u
    fl3 = fl2 ** u

    fl1 = frobenius_power(fl1, k, W, d, 1)
    fl2 = frobenius_power(fl2, k, W, d, 2)
    fl3 = frobenius_power(fl3, k, W, d, 3)
    fl4 = frobenius_power(fl4, k, W, d, 4)
    fl5 = frobenius_power(fl5, k, W, d, 5)
    fl6 = frobenius_power(fl6, k, W, d, 6)
    fl7 = frobenius_power(fl7, k, W, d, 7)

    M = fl5*fl7*fl0*fl2
    N = fl4*fl6*fl1*fl3
    N = frobenius_power(N, k, W, d, 12)

    t0 = M*N

    mk, sk, fk, sk_cyclo, ik = 162, 108, 22, 54, 343
    exp_u_n, mult, sq, inv, frob_power = 11, 13, 8, 1, 11

    if not NAF_rep:
        exp_u = (floor(log2(u)) - 1) * sk_cyclo + (hamming_weight(Integer(u).digits(2)) - 1) * mk
    else:
        exp_u = (ceil(log2(u)) - 1) * sk_cyclo + (NAf_hamming_weight(NAF(u)) - 1) * mk

    # count over Fp
    mult, sq, inv, frob = (mult * mk), (sq * sk), (inv * ik), (frob_power * fk)
    total_ops = (exp_u_n * exp_u) + mult + sq + inv + frob

    return t0, exp_u, mult, sq, inv, frob, total_ops
