def Line(D2, C):
    x2, y2, x22, x23 = D2[0], D2[1], D2[2], D2[3]
    l3, l2, l1, l0, l = C[0], C[1], C[2], C[3], C[5]
    cD2 = (y2 * l + l3 * x23 - l2 * x22 + l1 * x2 - l0)

    return cD2


def ADD(D1, D2, E, F):
    f0, f1, f2, f3, f4 = F[0], F[1], F[2], F[3], F[4]
    U11, U10, V11, V10 = D1[0], D1[1], D1[2], D1[3]
    U21, U20, V21, V20 = D2[0], D2[1], D2[2], D2[3]
    Z21, Z22, z21, z22 = D2[4], D2[5], D2[6], D2[7]
    z23 = Z21 * Z22
    z24 = z21 * z23
    Ut11, Ut10 = (U11 * z21), (U10 * z21)
    y1, y2 = (Ut11 - U21), (U20 - Ut10)
    y3 = U11 * y1
    y4 = y2 + y3
    r = y2 * y4 + y1 ^ 2 * U10
    inv1, inv0 = y1, y4
    w0 = V10 * z24 - V20
    w1 = V11 * z24 - V21
    w2 = inv0 * w0
    w3 = inv1 * w1
    s1p = y1 * w0 + y2 * w1
    s0p = w2 - U10 * w3
    rt = r * z23
    R = rt ^ 2
    Z31 = s1p * Z21
    Z32 = rt * Z21
    z31, z32 = (Z31 ^ 2), (Z32 ^ 2)
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
    U30 = s0tp * (s0p - 2 * s1p * U11) + s1p ^ 2 * (y3 - Ut10 - U20) + 2 * l1p + R * w1
    w1 = l2 * s1p - U31
    V30 = U30 * w1 - z31 * l0 * s1p
    V31 = U31 * w1 + z31 * (U30 - l1p)
    line = [l3, l2, l1, l0, l]
    D3 = [U31, U30, V31, V30, Z31, Z32, z31, z32]
    lcE = Line(E, line)

    return D3, lcE


def DBL(D1, E, F):
    f0, f1, f2, f3, f4 = F[0], F[1], F[2], F[3], F[4]
    U11, U10, V11, V10 = D1[0], D1[1], D1[2], D1[3]
    Z11, Z12, z11, z12 = D1[4], D1[5], D1[6], D1[7]
    Ut10 = U10 * z11
    Vt11 = 2 * V11
    Vt10 = 2 * V10
    zt11 = Vt10 * z11
    w0 = V11 ^ 2
    w1 = U11 ^ 2
    w2 = 4 * w0
    w3 = zt11 - U11 * Vt11
    r = Ut10 * w2 + zt11 * w3
    inv1p = -Vt11
    inv0p = w3
    z11p = z11 ^ 2
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
    R = rt ^ 2
    Z31 = s1p * Z11
    Z32 = rt * Z11
    z31 = Z31 ^ 2
    z32 = Z32 ^ 2
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
    lcE = Line(E, line)

    return D3, lcE


def new_coordinates(D):
    """
    :param D:
    :return:
    """
    u, v = D[0], D[1]
    # how to extract the coefficients in Sage?
    Cu, Cv = Coefficients(u), Coefficients(v)
    U1, U0 = Cu[1], Cu[0]
    V1, V0 = Cv[1], Cv[0]
    Z1, Z2 = (1), (1)
    z1, z2 = (Z1 ^ 2), (Z2 ^ 2)
    Dnew = [U1, U0, V1, V0, Z1, Z2, z1, z2]

    return Dnew


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
    # return 1
