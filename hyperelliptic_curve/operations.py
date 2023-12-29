def Line_case1(D2_vec, D2, C, L):
    """
    :param D2_vec:
    :param D2: [-xQ, yQ, xQ^2, -xQ^3]
    :param C: coefficients of the line
    :param L: 4-element vector
    :return: evaluation of the line
    """
    x2, y2, x22, x23 = D2[0], D2[1], D2_vec[0], D2_vec[1]
    l3, l2, l1, l0, l = C[0], C[1], C[2], C[3], C[4]
    c, c2, c3, c5 = L[0], L[1], L[2], L[3]

    l3, l2, l1, l0 = (l3 * c2), (l2 * c), (l1 * c3), (l0 * c5)
    cD2 = (y2 * l + l3 * x23 - l2 * x22 + l1 * x2 - l0)

    return cD2


def Line_case2(D2_vec, D2, D3, C):
    """
    :param D2_vec:
    :param D2:
    :param D3:
    :param C:
    :return:
    """
    # t1, t2, t3, t4, t5 = D2_vec[0], D2_vec[1], D2_vec[2], D2_vec[3], D2_vec[4]
    t6, t7, t8, t9, t10 = D2_vec[5], D2_vec[6], D2_vec[7], D2_vec[8], D2_vec[9]
    t11, t12, t13, t14, t15 = D2_vec[10], D2_vec[11], D2_vec[12], D2_vec[13], D2_vec[14]
    t16, t17, t18, t19, t20 = D2_vec[15], D2_vec[16], D2_vec[17], D2_vec[18], D2_vec[19]
    t21, t22, t23, t24, t25 = D2_vec[20], D2_vec[21], D2_vec[22], D2_vec[23], D2_vec[24]

    u21, u20 = D2[0][1], D2[0][0]
    # v21, v20 = D2[1][1], D2[1][0]

    # U31, U30, z31 = D3[0], D3[1], D3[6]
    l3, l2, l1, l0, l = C[0], C[1], C[2], C[3], C[4]

    w1, w2 = l, l3
    w3, w4 = (w1 * t6), (w2 * t17)
    w5, w6, w7 = (l2 * t12), (l1 * t9), (l0 * t8)
    w8 = w3 - w4 + w5 - w6 - w7
    w9 = w1 * w8
    w10, w11, w12, w13 = (w2 * t20), (l2 * t21), (l1 * t23), (l0 * t25)
    w14 = w10 - w11 + w12 - w13
    w15 = w2 * w14
    w16, w17, w18 = (l2 * t19), (l1 * t18), (l0 * t22)
    w19 = w16 - w17 + w18
    w20 = l2 * w19
    w21, w22 = (l1 * u20), (l0 * u21)
    w23 = w21 - w22
    w24, w25 = (l1 * w23), (l0 ** 2)

    cD2 = (w9 + w15 + w20 + w24 + w25)

    # i1, i2, i3 = (z31 ** 2), (z31 * U30), (U30 ** 2)
    # i4 = i1 * t19
    # i5, i6, i7 = (U31 * u20), (z31 * t18), (U30 * u21)
    # i8 = i5 - i6 - i7
    # i9, i10 = (U31 * i8), (i2 * t22)
    # u3D2 = i3 + i4 + i9 + i10

    lD2 = cD2

    return lD2


def precomputation_general_div(D2):
    """
    :param D2:
    :return:
    """
    u21, u20, v21, v20 = D2[0][1], D2[0][0], D2[1][1], D2[1][0]

    t1, t2 = u20 * v21, u21 * v20
    t3 = t1 - t2
    t4 = v21 * t3
    t5 = v20 ** 2
    t6 = t4 + t5
    t7 = u21 * v21
    t8 = 2 * v20 - t7
    t9, t10, t11 = (t1 + t3), (u21 * t3), (u20 * v20)
    t12, t13 = (t10 + 2 * t11), (u21 ** 2)
    t14 = t3 * t13
    t15 = 2 * t3 - t2
    t16 = u20 * t15
    t17 = t14 - t16
    t18, t19 = (u20 * u21), (u20 ** 2)
    t20 = t19 * u20
    t21 = t19 * u21
    t22 = t13 - 2 * u20
    t23 = u20 * t22
    t24 = t22 - u20
    t25 = u21 * t24

    Q = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24,
         t25]

    return Q


def ADD(D1, D2, Q_vec, Q, F, L=None, case: str = 'case1'):
    """
    :param D1:
    :param D2:
    :param Q_vec:
    :param Q:
    :param F:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :return:
    """
    # f0, f1, f2, f3, f4 = F[0], F[1], F[2], F[3], F[4]
    U11, U10, V11, V10 = D1[0], D1[1], D1[2], D1[3]
    U21, U20, V21, V20 = D2[0], D2[1], D2[2], D2[3]
    Z21, Z22, z21, z22 = D2[4], D2[5], D2[6], D2[7]
    z23 = Z21 * Z22
    z24 = z21 * z23
    Ut11, Ut10 = (U11 * z21), (U10 * z21)
    y1, y2 = (Ut11 - U21), (U20 - Ut10)
    y3 = U11 * y1
    y4 = y2 + y3

    r = y2 * y4 + y1 ** 2 * U10
    inv1, inv0 = y1, y4
    w0, w1 = (V10 * z24 - V20), (V11 * z24 - V21)
    w2, w3 = (inv0 * w0), (inv1 * w1)
    s1p = y1 * w0 + y2 * w1
    s0p = w2 - U10 * w3
    rt = r * z23
    R = rt ** 2
    Z31 = s1p * Z21
    Z32 = rt * Z21
    z31, z32 = (Z31 ** 2), (Z32 ** 2)
    s0tp = s0p * z21
    l2p = (s1p * U21)
    l2 = (l2p + s0tp)
    l0p = s0p * U20
    l0 = l0p + r * V20
    l1 = (s1p + s0p) * (U21 + U20) - l2p - l0p + r * V21
    l3, l = (s1p * z21), (rt * z21)
    w1 = Ut11 + U21
    U31 = s1p * (2 * s0tp - s1p * y1) - z32
    l1p = l1 * s1p
    U30 = s0tp * (s0p - 2 * s1p * U11) + s1p ** 2 * (y3 - Ut10 - U20) + 2 * l1p + R * w1
    w1 = l2 * s1p - U31
    V30 = U30 * w1 - z31 * l0 * s1p
    V31 = U31 * w1 + z31 * (U30 - l1p)
    line = [l3, l2, l1, l0, l]

    D3 = [U31, U30, V31, V30, Z31, Z32, z31, z32]
    if case == 'case1':
        lcE = Line_case1(Q_vec, Q, line, L)
    else:
        lcE = Line_case2(Q_vec, Q, D3, line)

    return D3, lcE


def DBL(D1, Q_vec, Q, F, L=None, case: str = 'case1'):
    """
    :param D1:
    :param Q_vec:
    :param Q:
    :param F:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :return:
    """
    f0, f1, f2, f3, f4 = F[0], F[1], F[2], F[3], F[4]
    U11, U10, V11, V10 = D1[0], D1[1], D1[2], D1[3]
    Z11, Z12, z11, z12 = D1[4], D1[5], D1[6], D1[7]
    Ut10 = U10 * z11
    Vt11 = 2 * V11
    Vt10 = 2 * V10
    zt11 = Vt10 * z11
    w0 = V11 ** 2
    w1 = U11 ** 2
    w2 = 4 * w0
    w3 = zt11 - U11 * Vt11
    r = Ut10 * w2 + zt11 * w3
    inv1p = -Vt11
    inv0p = w3
    z11p = z11 ** 2
    z11pp = z11 * z11p
    w3 = f3 * z11p + w1
    w4 = 2 * Ut10
    k1p = z12 * (2 * w1 + w3 - w4)
    k0p = z12 * (U11 * (2 * w4 - w3) + f2 * z11pp) - w0
    w0 = k0p * inv0p
    w1 = k1p * inv1p
    s1p = z11 * (zt11 * k1p - Vt11 * k0p)
    s0p = w0 - Ut10 * w1
    z13 = Z11 * Z12
    rt = r * z13
    R = rt ** 2
    Z31 = s1p * Z11
    Z32 = rt * Z11
    z31 = Z31 ** 2
    z32 = Z32 ** 2
    st0p = s0p * z11
    l2p = s1p * U11
    l2 = l2p + st0p
    l0p = s0p * U10
    l0 = l0p + r * V10
    rp = r * V11
    l1 = (s1p + s0p) * (U11 + U10) - l2p - l0p + rp
    l3 = s1p * z11
    l = rt * z11
    U30 = 2 * (rp * s1p + R * U11) + s0p * st0p
    U31 = 2 * s1p * st0p - z32
    w1 = l2 * s1p - U31
    V30 = U30 * w1 - z31 * l0 * s1p
    V31 = U31 * w1 + z31 * (U30 - l1 * s1p)
    line = [l3, l2, l1, l0, l]

    D3 = [U31, U30, V31, V30, Z31, Z32, z31, z32]
    if case == 'case1':
        lcE = Line_case1(Q_vec, Q, line, L)
    else:
        lcE = Line_case2(Q_vec, Q, D3, line)

    return D3, lcE


def new_coordinates(D) -> list:
    """
    :param D: [x^2+U1x+U0, V1x + V0]
    :return: 8-element vector
    """
    U1, U0 = D[0][1], D[0][0]
    V1, V0 = D[1][1], D[1][0]
    Z1, Z2 = (1), (1)
    z1, z2 = (Z1 ** 2), (Z2 ** 2)
    D_ = [U1, U0, V1, V0, Z1, Z2, z1, z2]

    return D_


def HEC_random_point(C):
    """
    :param C:
    :return:
    """
    f = C.hyperelliptic_polynomials()[0]
    while True:
        x_r = f.base_ring().random_element()
        y2 = f(x_r)
        if y2.is_square():
            return [x_r, y2.sqrt()]


def HEC_random_points_uniq(C, n):
    """
    :param C:
    :param n:
    :return:
    """
    f = C.hyperelliptic_polynomials()[0]
    res = []
    for i in range(1, n + 1):
        tries = 100
        found = False
        while tries > 0 and (not found):
            P = HEC_random_point(C)
            if not (P in res):
                res.append(P)
                found = True
            tries = tries - 1

    return res


def JC_random_element(C):
    """
    :param C:
    :return:
    """
    f = C.hyperelliptic_polynomials()[0]
    R = f.base_ring()['x']
    (x,) = R._first_ngens(1)
    J = C.jacobian()
    points = HEC_random_points_uniq(C, C.genus())
    u = 1
    for point in points:
        u = u * (x - point[0])
    v = f.parent().lagrange_polynomial(points)

    return J([u, v])
