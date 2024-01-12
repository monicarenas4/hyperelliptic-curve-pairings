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
    f = f.frobenius(4) / f  # ---> 1 f + 1 m + 1 i : 1 frobenius + 1 inversion + 1 multiplication
    fp = f.frobenius(1)  # ---> 1 f
    fu1 = f ** u  # ---> 1 e_u : 1 exponentiation by u
    f1 = fp * fu1  # ---> 1 m
    fp2 = f1.frobenius(2)  # ---> 1 f
    f1u1 = f1 ** u  # ---> 1 e_u
    f1u2 = f1u1 ** u  # ---> 1 e_u
    f2 = fp2 * f1u2  # ---> 1 m
    f22 = f2 ** 2  # ---> 1 s
    y1 = f22 ** u0  # ---> 1 e_u0
    y2 = y1 ** u0  # ---> 1 e_u0
    y3 = y2 ** u0  # ---> 1 e_u0
    y4 = y3 ** u0  # ---> 1 e_u0
    y32 = y3 ** 2  # ---> 1 s
    y34 = y32 ** 2  # ---> 1 s
    y38 = y34 ** 2  # ---> 1 s
    y42 = y4 ** 2  # ---> 1 s
    y44 = y42 ** 2  # ---> 1 s
    y48 = y44 ** 2  # ---> 1 s
    y416 = y48 ** 2  # ---> 1 s
    y4162 = y416 ** 2  # ---> 1 s
    y4164 = y4162 ** 2  # ---> 1 s
    y4168 = y4164 ** 2  # ---> 1 s
    y41616 = y4168 ** 2  # ---> 1 s
    y41624 = y41616 * y4168  # ---> 1 m
    y41625 = y41624 * y416  # ---> 1 m
    y382 = y38 ** 2  # ---> 1 s
    y384 = y382 ** 2  # ---> 1 s
    y385 = y384 * y38  # ---> 1 m
    x1 = y4168 * f2  # ---> 1 m
    x2 = x1 ** lx  # 1 e_lx : 1 exponentiation by lx
    x3 = x2 * y384  # ---> 1 m
    x4 = x3 ** lx  # ---> 1 e_lx : 1 exponentiation by lx
    z1 = y41625 * y385  # ---> 1 m
    z2 = z1 * f22  # ---> 1 m
    N = f * x4 * z2  # ---> 2 m
    y22 = y2 ** 2  # ---> 1 s
    y23 = y22 * y2  # ---> 1 m
    y22lx = y22 ** lx  # ---> 1 e_lx : 1 exponentiation by lx
    m = y22lx * y2 * y1  # ---> 2 m
    m2 = m ** 2  # ---> 1 s
    m4 = m2 ** 2  # ---> 1 s
    M = m4 * y23  # ---> 1 m
    t0 = N / M  # ---> 1 m + 1 i : 1 inversion + 1 multiplication

    exp_u, exp_u0, exp_lx, mult, sq, inv, frob = 3, 4, 3, 17, 17, 2, 3

    return t0, exp_u, exp_u0, exp_lx, mult, sq, inv, frob


def final_exponentiation_k16(f, U, W, k=16):
    """
    :param f:
    :param U:
    :param W:
    :param k:
    :return:
    """
    u, um = U[0], U[1]
    f = frobenius_power(f, k, W, 8) / f  # ---> 1 f + 1 m + 1 i
    f1 = f
    f2 = f1 ** 2  # ---> 1 s
    f4 = f2 ** 2  # ---> 1 s
    f8 = f4 ** 2  # ---> 1 s
    fum1 = f ** um  # ? ---> 1 e_um : 1 exponentiation by um (I will change the names of the variables at some point)
    fu1 = fum1 * f1  # ---> 1 m
    # fup1 = fu1 * f1  # fup1 not used
    fum2 = fum1 ** um  # ---> 1 e_um : 1 exponentiation by um
    f2u1 = fu1 ** 2  # ---> 1 s
    f4u1 = f2u1 ** 2  # ---> 1 s
    fup2 = fum2 * f4u1  # ---> 1 m
    g1 = fum2 ** u  # ---> 1 e_u
    g2 = g1 ** u  # ---> 1 e_u
    g3 = g2 ** u  # ---> 1 e_u
    g4 = g3 ** u  # ---> 1 e_u
    g5 = g4 ** 2  # ---> 1 s
    fl0 = fup2 * g5  # ---> 1 m
    fl0u = fl0 ** u  # ---> 1 e_u
    fl0u2 = fl0u ** u  # ---> 1 e_u
    fl0u3 = fl0u2 ** u  # ---> 1 e_u
    fn = f8 * fl0u3  # ---> 1 m
    fnu1 = fn ** u  # ---> 1 e_u
    fnu2 = fnu1 ** u  # ---> 1 e_u
    fnu3 = fnu2 ** u  # ---> 1 e_u
    fnu4 = fnu3 ** u  # # ---> 1 e_u
    N0 = fl0
    N1 = frobenius_power(fn, k, W, 1)  # ---> 1 f
    N2 = frobenius_power(fnu3, k, W, 2)  # ---> 1 f
    N3 = frobenius_power(fl0u, k, W, 3)  # ---> 1 f
    N4 = frobenius_power(fnu1, k, W, 4)  # ---> 1 f
    N5 = frobenius_power(fnu4, k, W, 5)  # ---> 1 f
    N6 = frobenius_power(fl0u2, k, W, 6)  # ---> 1 f
    N7 = frobenius_power(fnu2, k, W, 7)  # ---> 1 f

    N, M = (N0 * N2 * N3 * N5 * N6), (N1 * N4 * N7)  # ---> 6 m
    t0 = N / M  # ---> 1 m + 1 i

    exp_u, exp_um, mult, sq, inv, frob_power = 2, 11, 12, 6, 2, 8

    return t0, exp_u, exp_um, mult, sq, inv, frob_power
