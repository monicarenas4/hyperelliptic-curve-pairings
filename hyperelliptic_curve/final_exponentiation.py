def final_exponentiation(f, U, K):
    """
    :param f:
    :param U:
    :param K: field Fp
    :return:
    """
    u, u0, lx, ly = U[0], U[1], U[2], U[3]
    t0 = 1
    # f = Frobenius(f, K, 4) / f
    f = f.frobenius(4) / f
    # fp = Frobenius(f, K, 1)
    fp = f.frobenius(1)
    fu1 = f ^ u
    f1 = fp * fu1
    # fp2 = Frobenius(f1, K, 2)
    fp2 = f1.frobenius(2)
    f1u1 = f1 ^ u
    f1u2 = f1u1 ^ u
    f2 = fp2 * f1u2
    f22 = f2 ^ 2
    y1 = f22 ^ u0
    y2 = y1 ^ u0
    y3 = y2 ^ u0
    y4 = y3 ^ u0
    y32 = y3 ^ 2
    y34 = y32 ^ 2
    y38 = y34 ^ 2
    y42 = y4 ^ 2
    y44 = y42 ^ 2
    y48 = y44 ^ 2
    y416 = y48 ^ 2
    y4162 = y416 ^ 2
    y4164 = y4162 ^ 2
    y4168 = y4164 ^ 2
    y41616 = y4168 ^ 2
    y41624 = y41616 * y4168
    y41625 = y41624 * y416
    y382 = y38 ^ 2
    y384 = y382 ^ 2
    y385 = y384 * y38
    x1 = y4168 * f2
    x2 = x1 ^ lx
    x3 = x2 * y384
    x4 = x3 ^ lx
    z1 = y41625 * y385
    z2 = z1 * f22
    N = f * x4 * z2
    y22 = y2 ^ 2
    y23 = y22 * y2
    y22lx = y22 ^ lx
    m = y22lx * y2 * y1
    m2 = m ^ 2
    m4 = m2 ^ 2
    M = m4 * y23
    t0 = N / M

    return t0

