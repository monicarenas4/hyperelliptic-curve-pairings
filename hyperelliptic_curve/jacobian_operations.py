def line_case1(D2_vec, D2, C, L, twist: str = None):
    """
    :param D2_vec:
    :param D2: [-xQ, yQ, xQ^2, -xQ^3]
    :param C: coefficients of the line
    :param L: 4-element vector
    :return: evaluation of the line
    """

    l3, l2, l1, l0, l = C[0], C[1], C[2], C[3], C[4]
    c, c2, c3, c5 = L[0], L[1], L[2], L[3]

    if twist == None:
        x2, y2, x22, x23 = -D2[0], D2[1], D2_vec[0], D2_vec[1]
        x2 = x2 / c ** 2
        y2 = y2 / c ** 5
        x22 = x22 / c ** 4
        x23 = x23 / c ** 6
        l3, l2, l1, l0 = l3, l2, l1, l0
    else:
        x2, y2, x22, x23 = -D2[0], D2[1], D2_vec[0], D2_vec[1]
        l3, l2, l1, l0 = (l3 * c), (l2 / c), (l1 / c ** 3), (l0 / c ** 5)

    cD2 = (y2 * l + l3 * x23 - l2 * x22 + l1 * x2 - l0)  # M4

    return cD2


def line_case2(D2_vec, D2, D3, C, L, twist=None):
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

    c, c2, c3, c5 = L[0], L[1], L[2], L[3]

    u21, u20 = D2[0][1], D2[0][0]
    # v21, v20 = D2[1][1], D2[1][0]

    # U31, U30, z31 = D3[0], D3[1], D3[6]
    l3, l2, l1, l0, l = C[0], C[1], C[2], C[3], C[4]

    w1, w2 = l, l3
    w3, w4 = (w1 * t6), (w2 * t17)  # M2
    w5, w6, w7 = (l2 * t12), (l1 * t9), (l0 * t8)  # M5

    if twist == None:
        w8 = w3 / (c ** 10) - w4 / (c ** 11) + w5 / (c ** 9) - w6 / (c ** 7) - w7 / (c ** 5)
        w9 = w1 * w8  # M6
        w10, w11, w12, w13 = (w2 * t20), (l2 * t21), (l1 * t23), (l0 * t25)  # M10
        w14 = w10 / (c ** 12) - w11 / (c ** 10) + w12 / (c ** 8) - w13 / (c ** 6)
        w15 = w2 * w14  # M11
        w16, w17, w18 = (l2 * t19), (l1 * t18), (l0 * t22)  # M14
        w19 = w16 / (c ** 8) - w17 / (c ** 6) + w18 / (c ** 4)
        w20 = l2 * w19  # M15
        w21, w22 = (l1 * u20), (l0 * u21)  # M17
        w23 = w21 / (c ** 4) - w22 / (c ** 2)
        w24, w25 = (l1 * w23), (l0 ** 2)  # M18 S1
        cD2 = (w9 + w15 + w20 + w24 + w25)
    else:
        w8 = w3 / (c ** 2) - w4 / c + w5 / (c ** 3) - w6 / (c ** 5) - w7 / (c ** 7)
        w9 = w1 * w8
        w10, w11, w12, w13 = (w2 * t20), (l2 * t21), (l1 * t23), (l0 * t25)
        w14 = w10 - w11 / (c ** 2) + w12 / (c ** 4) - w13 / (c ** 6)
        w15 = w2 * w14
        w16, w17, w18 = (l2 * t19), (l1 * t18), (l0 * t22)
        w19 = w16 / (c ** 4) - w17 / (c ** 6) + w18 / (c ** 8)
        w20 = l2 * w19
        w21, w22 = (l1 * u20), (l0 * u21)
        w23 = w21 / (c ** 8) - w22 / (c ** 10)
        w24, w25 = (l1 * w23), (l0 ** 2)
        cD2 = (w9 + w15 + w20 + w24 + w25 / (c ** 12))

    # i1, i2, i3 = (z31 ** 2), (z31 * U30), (U30 ** 2)
    # i4 = i1 * t19
    # i5, i6, i7 = (U31 * u20), (z31 * t18), (U30 * u21)
    # i8 = i5 - i6 - i7
    # i9, i10 = (U31 * i8), (i2 * t22)
    # u3D2 = i3 + i4 + i9 + i10
    # lD2 = cD2

    return cD2

def precomputation_degenerate_div(P):
    """
    :param P:
    :return:
    """
    xP, yP = P[0], P[1]
    a = xP ** 2  # S1
    b = -a * xP  # M1
    P_prec = [a, b]

    return P_prec


def precomputation_general_div(D2):
    """
    :param D2:
    :return:
    """
    u21, u20, v21, v20 = D2[0][1], D2[0][0], D2[1][1], D2[1][0]

    t1 = u20 * v21  # M1
    t2 = u21 * v20  # M2
    t3 = t1 - t2
    t4 = v21 * t3  # M3
    t5 = v20 ** 2  # S1
    t6 = t4 + t5
    t7 = u21 * v21  # M4
    t8 = 2 * v20 - t7
    t9 = t1 + t3
    t10 = u21 * t3  # M5
    t11 = u20 * v20  # M6
    t12 = t10 + 2 * t11
    t13 = u21 ** 2  # S2
    t14 = t3 * t13  # M7
    t15 = 2 * t3 - t2
    t16 = u20 * t15  # M8
    t17 = t14 - t16
    t18 = u20 * u21  # M9
    t19 = u20 ** 2  # S3
    t20 = t19 * u20  # M10
    t21 = t19 * u21  # M11
    t22 = t13 - 2 * u20
    t23 = u20 * t22  # M12
    t24 = t22 - u20
    t25 = u21 * t24  # M13

    Q = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24,
         t25]

    return Q


def ADD(D1, D2, Q_vec, Q, F, L=None, case: str = 'case1', twist: str = None):
    """
    :param D1:
    :param D2:
    :param Q_vec:
    :param Q:
    :param F:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :return:
    """
    f0, f1, f2, f3, f4 = F[0], F[1], F[2], F[3], F[4]
    U11, U10, V11, V10 = D1[0], D1[1], D1[2], D1[3]
    U21, U20, V21, V20 = D2[0], D2[1], D2[2], D2[3]
    Z21, Z22, z21, z22 = D2[4], D2[5], D2[6], D2[7]
    z23 = Z21 * Z22  # M1
    z24 = z21 * z23  # M2
    Ut11, Ut10 = (U11 * z21), (U10 * z21)  # M4
    y1, y2 = (Ut11 - U21), (U20 - Ut10)
    y3 = U11 * y1  # M5
    y4 = y2 + y3

    r = y2 * y4 + y1 ** 2 * U10  # M6 S1
    inv1, inv0 = y1, y4
    w0, w1 = (V10 * z24 - V20), (V11 * z24 - V21)  # M8
    w2, w3 = (inv0 * w0), (inv1 * w1)  # M10
    s1p = y1 * w0 + y2 * w1  # M12
    s0p = w2 - U10 * w3  # M13
    rt = r * z23  # M14
    R = rt ** 2  # S2
    Z31 = s1p * Z21  # M15
    Z32 = rt * Z21  # M16
    z31, z32 = (Z31 ** 2), (Z32 ** 2)  # S4
    s0tp = s0p * z21  # M17
    l2p = (s1p * U21)  # M18
    l2 = (l2p + s0tp)
    l0p = s0p * U20  # M19
    l0 = l0p + r * V20  # M20
    l1 = (s1p + s0p) * (U21 + U20) - l2p - l0p + r * V21  # M22
    l3, l = (s1p * z21), (rt * z21)  # M24
    w1 = Ut11 + U21
    U31 = s1p * (2 * s0tp - s1p * y1) - z32  # M26
    l1p = l1 * s1p  # M27
    U30 = s0tp * (s0p - 2 * s1p * U11) + s1p ** 2 * (y3 - Ut10 - U20) + 2 * l1p + R * w1  # M30 S5
    w1 = l2 * s1p - U31  # M31
    V30 = U30 * w1 - z31 * l0 * s1p  # M34
    V31 = U31 * w1 + z31 * (U30 - l1p)  # M36
    line = [l3, l2, l1, l0, l]

    D3 = [U31, U30, V31, V30, Z31, Z32, z31, z32]
    if case == 'case1' and twist != None:
        lcE = line_case1(Q_vec, Q, line, L, twist='k16')
    elif case == 'case1':
        lcE = line_case1(Q_vec, Q, line, L)
    elif case == 'case2' and twist != None:
        lcE = line_case2(Q_vec, Q, D3, line, L, twist='k16')
    else:
        lcE = line_case2(Q_vec, Q, D3, line, L)

    return D3, lcE


def DBL(D1, Q_vec, Q, F, L=None, case: str = 'case1', twist: str = None):
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
    Ut10 = U10 * z11  # M1
    Vt11 = 2 * V11
    Vt10 = 2 * V10
    zt11 = Vt10 * z11  # M2
    w0 = V11 ** 2  # S1
    w1 = U11 ** 2  # S2
    w2 = 4 * w0
    w3 = zt11 - U11 * Vt11  # M3
    r = Ut10 * w2 + zt11 * w3  # M5
    inv1p = -Vt11
    inv0p = w3
    z11p = z11 ** 2  # S3
    z11pp = z11 * z11p  # M6
    w3 = f3 * z11p + w1  # f3=0
    w4 = 2 * Ut10
    k1p = z12 * (2 * w1 + w3 - w4)  # M7
    k0p = z12 * (U11 * (2 * w4 - w3) + f2 * z11pp) - w0  # M9 f2 =0
    w0 = k0p * inv0p  # M10
    w1 = k1p * inv1p  # M11
    s1p = z11 * (zt11 * k1p - Vt11 * k0p)  # M14
    s0p = w0 - Ut10 * w1  # M15
    z13 = Z11 * Z12  # M16
    rt = r * z13  # M17
    R = rt ** 2  # S4
    Z31 = s1p * Z11  # M18
    Z32 = rt * Z11  # M19
    z31 = Z31 ** 2  # S5
    z32 = Z32 ** 2  # S6
    st0p = s0p * z11  # M20
    l2p = s1p * U11  # M21
    l2 = l2p + st0p
    l0p = s0p * U10  # M22
    l0 = l0p + r * V10  # M23
    rp = r * V11  # M24
    l1 = (s1p + s0p) * (U11 + U10) - l2p - l0p + rp  # M25
    l3 = s1p * z11  # M26
    l = rt * z11  # M27
    U30 = 2 * (rp * s1p + R * U11) + s0p * st0p  # M30
    U31 = 2 * s1p * st0p - z32  # M31
    w1 = l2 * s1p - U31  # M32
    V30 = U30 * w1 - z31 * l0 * s1p  # M35
    V31 = U31 * w1 + z31 * (U30 - l1 * s1p)  # M38
    line = [l3, l2, l1, l0, l]

    D3 = [U31, U30, V31, V30, Z31, Z32, z31, z32]

    if case == 'case1' and twist != None:
        lcE = line_case1(Q_vec, Q, line, L, twist='k16')
    elif case == 'case1':
        lcE = line_case1(Q_vec, Q, line, L)
    elif case == 'case2' and twist != None:
        lcE = line_case2(Q_vec, Q, D3, line, L, twist='k16')
    else:
        lcE = line_case2(Q_vec, Q, D3, line, L)

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
    f = C.hyperelliptic_polynomials()[0]
    while True:
        x_r = f.base_ring().random_element()
        y2 = f(x_r)
        if y2.is_square():
            return [x_r, y2.sqrt()]


def HEC_random_points_uniq(C, n):
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
