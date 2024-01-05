from _utils import frobenius_power


def final_exponentiation_cp8(f, U: list, K):
    """
    :param f: output of the miller loop
    :param U: 4-elements vector
    :param K:
    :return: final pairing output
    """
    u, u0, lx, ly = U[0], U[1], U[2], U[3]
    t0 = 1
    f = f.frobenius(4) / f
    fp = f.frobenius(1)
    fu1 = f ** u
    f1 = fp * fu1
    # fp2 = Frobenius(f1, K, 2)
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
    N = f * x4 * z2
    y22 = y2 ** 2
    y23 = y22 * y2
    y22lx = y22 ** lx
    m = y22lx * y2 * y1
    m2 = m ** 2
    m4 = m2 ** 2
    M = m4 * y23
    t0 = N / M

    return t0


def final_exponentiation_k16(f, U, W, k=16):
    """
    :param f:
    :param U:
    :param W:
    :param k:
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
    fup1 = fu1 * f1
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

    return t0
