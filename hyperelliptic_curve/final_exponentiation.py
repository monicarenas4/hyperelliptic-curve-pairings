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
    fu1 = f ** u  # S1
    f1 = fp * fu1  # M1
    # fp2 = Frobenius(f1, K, 2)
    fp2 = f1.frobenius(2)
    f1u1 = f1 ** u  # S2
    f1u2 = f1u1 ** u  # S3
    f2 = fp2 * f1u2  # M2
    f22 = f2 ** 2  # S4
    y1 = f22 ** u0  # S5
    y2 = y1 ** u0  # S6
    y3 = y2 ** u0  # S7
    y4 = y3 ** u0  # S8
    y32 = y3 ** 2  # S9
    y34 = y32 ** 2  # S10
    y38 = y34 ** 2  # S11
    y42 = y4 ** 2  # S12
    y44 = y42 ** 2  # S13
    y48 = y44 ** 2  # S14
    y416 = y48 ** 2  # S15
    y4162 = y416 ** 2  # S16
    y4164 = y4162 ** 2  # S17
    y4168 = y4164 ** 2  # S18
    y41616 = y4168 ** 2  # S19
    y41624 = y41616 * y4168  # M3
    y41625 = y41624 * y416  # M4
    y382 = y38 ** 2  # S20
    y384 = y382 ** 2  # S21
    y385 = y384 * y38  # M5
    x1 = y4168 * f2  # M6
    x2 = x1 ** lx  # S22
    x3 = x2 * y384  # M7
    x4 = x3 ** lx  # S23
    z1 = y41625 * y385  # M8
    z2 = z1 * f22  # M9
    N = f * x4 * z2  # M11
    y22 = y2 ** 2  # S24
    y23 = y22 * y2  # M12
    y22lx = y22 ** lx  # S25
    m = y22lx * y2 * y1  # M14
    m2 = m ** 2  # S26
    m4 = m2 ** 2  # S27
    M = m4 * y23  # M15
    t0 = N / M

    return t0, 15, 27


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
    fnu4 = fnu3 ** u  # S19
    N0 = fl0
    N1 = frobenius_power(fn, k, W, 1)
    N2 = frobenius_power(fnu3, k, W, 2)
    N3 = frobenius_power(fl0u, k, W, 3)
    N4 = frobenius_power(fnu1, k, W, 4)
    N5 = frobenius_power(fnu4, k, W, 5)
    N6 = frobenius_power(fl0u2, k, W, 6)
    N7 = frobenius_power(fnu2, k, W, 7)

    N, M = (N0 * N2 * N3 * N5 * N6), (N1 * N4 * N7)  # M11
    t0 = N / M

    return t0, 11, 19
