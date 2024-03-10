def line_case1(Q_prec, Q, line, c_vec, twist: bool = False):
    """
    :param Q_prec: precomputation vector needed for line evaluation
    :param Q: a point on the curve Ct(Fp^s); line will be evaluated at Q
    :param line: coefficients of the line to be evaluated
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param twist: distinguish between twisted ate and ate_i pairings
    :return: evaluation of the line c(Q)

    Evaluate line in case where 2nd pairing input is a degenerate divisor
    See formulas in Fan et al. https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf

    cost: 4 mult.
    """

    # line coefficients
    # line 1 refers to case 1, i.e. the 2nd pairing input is a degenerate divisor
    l3, l2, l1, l0, l = line[0], line[1], line[2], line[3], line[4]

    if not twist:
        x2, y2, x22, x23 = -Q[0], Q[1], Q_prec[0], Q_prec[1]
        x2 = x2 * c_vec[2]
        y2 = y2 * c_vec[5]
        x22 = x22 * c_vec[4]
        x23 = x23 * c_vec[6]
        l3, l2, l1, l0 = l3, l2, l1, l0
    else:
        x2, y2, x22, x23 = -Q[0], Q[1], Q_prec[0], Q_prec[1]
        # transition from Fq = Fp^s to Fp^k = Fq^8
        l3, l2, l1, l0 = (l3 * c_vec[0]), (l2 * c_vec[1]), (l1 * c_vec[3]), (l0 * c_vec[5])

    # line evaluation
    cQ = (y2 * l + l3 * x23 - l2 * x22 + l1 * x2 - l0)  # M4
    const_mult_, mult_, sq_ = 4, 0, 0

    return cQ, const_mult_, mult_, sq_


def line_case2(Q_prec, Q, line, c_vec, twist: bool = False):
    """
    :param Q_prec: precomputation vector needed for line evaluation
    :param Q: a general divisor in Jt(Fp^s); line will be evaluated at Q
    :param D3: not used because of the twist
    :param line: coefficients of the line to be evaluated
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param twist: distinguish between twisted ate and ate_i pairings
    :return: evaluation of the line c(Q)

    Evaluate line in case where 2nd pairing input is a general divisor
    See formulas in Fan et al. https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf

    cost: 14 constant mult. + 14 mult. + 1 sq.
    """

    # rename precomputation values to agree with Fan et al. formulas
    t6, t8, t9, t12 = Q_prec[0], Q_prec[1], Q_prec[2], Q_prec[3]
    t17, t18, t19, t20 = Q_prec[4], Q_prec[5], Q_prec[6], Q_prec[7]
    t21, t22, t23, t25 = Q_prec[8], Q_prec[9], Q_prec[10], Q_prec[11]

    u21, u20 = Q[0][1], Q[0][0]

    # line coefficients
    # line 1 refers to case 1, i.e. the 2nd pairing input is a degenerate divisor
    l3, l2, l1, l0, l = line[0], line[1], line[2], line[3], line[4]

    # line evaluation
    w1, w2 = l, l3
    w3, w4 = (w1 * t6), (w2 * t17)
    w5, w6, w7 = (l2 * t12), (l1 * t9), (l0 * t8)

    if not twist:
        w8 = w3 * c_vec[10] - w4 * c_vec[11] + w5 * c_vec[9] - w6 * c_vec[7] - w7 * c_vec[5] # transition from Fq = Fp^s to Fp^k = Fq^8
        w9 = w1 * w8
        w10, w11, w12, w13 = (w2 * t20), (l2 * t21), (l1 * t23), (l0 * t25)
        w14 = w10 * c_vec[12] - w11 * c_vec[10] + w12 * c_vec[8] - w13 * c_vec[6] # transition from Fq = Fp^s to Fp^k = Fq^8
        w15 = w2 * w14
        w16, w17, w18 = (l2 * t19), (l1 * t18), (l0 * t22)
        w19 = w16 * c_vec[8] - w17 * c_vec[6] + w18 * c_vec[4] # transition from Fq = Fp^s to Fp^k = Fq^8
        w20 = l2 * w19
        w21, w22 = (l1 * u20), (l0 * u21)
        w23 = w21 * c_vec[4] - w22 * c_vec[2]  # transition from Fq = Fp^s to Fp^k = Fq^8
        w24, w25 = (l1 * w23), (l0 ** 2)
        # line evaluation
        cQ = (w9 + w15 + w20 + w24 + w25)
    else:
        w8 = w3 * c_vec[2] - w4 * c_vec[1] + w5 * c_vec[3] - w6 * c_vec[5] - w7 * c_vec[7] # transition from Fq = Fp^s to Fp^k = Fq^8
        w9 = w1 * w8
        w10, w11, w12, w13 = (w2 * t20), (l2 * t21), (l1 * t23), (l0 * t25)
        w14 = w10 - w11 * c_vec[2] + w12 * c_vec[4] - w13 * c_vec[6] # transition from Fq = Fp^s to Fp^k = Fq^8
        w15 = w2 * w14
        w16, w17, w18 = (l2 * t19), (l1 * t18), (l0 * t22)
        w19 = w16 * c_vec[4] - w17 * c_vec[6] + w18 * c_vec[8] # transition from Fq = Fp^s to Fp^k = Fq^8
        w20 = l2 * w19
        w21, w22 = (l1 * u20), (l0 * u21)
        w23 = w21 * c_vec[8] - w22 * c_vec[10] # transition from Fq = Fp^s to Fp^k = Fq^8
        w24, w25 = (l1 * w23), (l0 ** 2)
        # line evaluation
        cQ = (w9 + w15 + w20 + w24 + w25 * c_vec[12]) # transition from Fq = Fp^s to Fp^k = Fq^8

    const_mult_, mult_, sq_ = 14, 14, 1

    return cQ, const_mult_, mult_, sq_


def precomputation_degenerate_div(P):
    """
    :param P: a point on the curve C s.t. P = (xP, yP)
    :return: precomputation of values needed for line 1 evaluation: P_prec = [xP^2, - xP^3]

    cost: 1 mult. + 1 sq.
    """
    # Precomputation of values needed for line evaluation in case where 2nd pairing input is a degenerate divisor
    # See formulas in Fan et al. https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf
    xP, yP = P[0], P[1]
    a = xP ** 2
    b = -a * xP

    # precomputation vector (2 elements)
    P_prec = [a, b]
    mult_, sq_ = 1, 1

    return P_prec, mult_, sq_

def precomputation_general_div(P):
    """
    :param P: a general divisor in the Jacobian J s.t. P = [x^2 + u21*x + u21, v21*x + v20]
    :return:

    cost: 13 mult. + 3 sq.
    """
    # Precomputation of values needed for line evaluation in case where 2nd pairing input is a general divisor
    # See formulas in Fan et al. https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf
    u21, u20, v21, v20 = P[0][1], P[0][0], P[1][1], P[1][0]

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

    # precomputation vector (25 elements)
    P_prec = [t6, t8, t9, t12, t17, t18, t19, t20, t21, t22, t23, t25]

    mult_, sq_ = 13, 3

    return P_prec, mult_, sq_


def ADD(P, T, Q_prec, Q, F, c_vec, case: str = 'case1', twist: bool = False):
    """
    :param P: general divisor P = [x^2 + U11*x + U10, V11*x + V10]
    :param T: divisor in Fan et al. coordinate system T = [U21, U20, V21, V20, Z21, Z22, z21, z22]
    :param Q_prec: precomputation vector related to divisor Q
    :param Q: degenerate or general divisor for line evaluation
    :param F: curve coefficients C/Fp: y^2 = f(x)
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param twist: twist = True => P in Jt(Fp^s), Q in J(Fp) or twist = False => P in J(Fp), Q in Jt(Fp^s)
    :return: divisor: R = T + P (in Fan et al. coordinate system) and lcE: line evaluated at Q

    Add two divisors T, P, where T = [U21, U20, V21, V20, Z21, Z22, z21, z22] and P = [x^2 + U11*x + U10, V11*x + V10]
    Evaluate line at a degenerate or general divisor Q
    See formulas in Fan et al. https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf

    cost: 37 mult. + 5 sq.
    """

    U11, U10, V11, V10 = P[0], P[1], P[2], P[3]
    U21, U20, V21, V20 = T[0], T[1], T[2], T[3]
    Z21, Z22, z21, z22 = T[4], T[5], T[6], T[7]
    z23 = Z21 * Z22  # 1m
    z24 = z21 * z23  # 1m
    Ut11, Ut10 = (U11 * z21), (U10 * z21)  # 2m
    y1, y2 = (Ut11 - U21), (U20 - Ut10)
    y3 = U11 * y1  # 1m
    y4 = y2 + y3

    r = y2 * y4 + y1 ** 2 * U10  # 2m, 1s
    inv1, inv0 = y1, y4
    w0, w1 = (V10 * z24 - V20), (V11 * z24 - V21)  # M8 2m
    w2, w3 = (inv0 * w0), (inv1 * w1)  # M10 2m
    s1p = y1 * w0 + y2 * w1  # M12 2m
    s0p = w2 - U10 * w3  # M13 1m
    rt = r * z23  # M14 1m
    R = rt ** 2  # S2
    Z31 = s1p * Z21  # M15 1m
    Z32 = rt * Z21  # M16 1m
    z31, z32 = (Z31 ** 2), (Z32 ** 2)  # S4
    s0tp = s0p * z21  # M17 1m
    l2p = (s1p * U21)  # M18 1m
    l2 = (l2p + s0tp)
    l0p = s0p * U20  # M19 1m
    l0 = l0p + r * V20  # M20 1m
    l1 = (s1p + s0p) * (U21 + U20) - l2p - l0p + r * V21  # M22 2m
    l3, l = (s1p * z21), (rt * z21)  # M24 2m
    w1 = Ut11 + U21
    U31 = s1p * (2 * s0tp - s1p * y1) - z32  # M26 2m
    l1p = l1 * s1p  # M27 1m
    U30 = s0tp * (s0p - 2 * s1p * U11) + s1p ** 2 * (y3 - Ut10 - U20) + 2 * l1p + R * w1  # M31 S5 4m
    w1 = l2 * s1p - U31  # M32 1m
    V30 = U30 * w1 - z31 * l0 * s1p  # M35 3m
    V31 = U31 * w1 + z31 * (U30 - l1p)  # M37 2m
    line = [l3, l2, l1, l0, l]

    # cost for addition + line computation: 37 mult. + 5 sq.
    mult_, sq_ = 37, 5

    R = [U31, U30, V31, V30, Z31, Z32, z31, z32]
    if case == 'case1' and twist is True:
        # line 1: for degenerate divisor Q in Jt(Fp^s)
        lcE, const_mult_line, mult_line, sq_line = line_case1(Q_prec, Q, line, c_vec, twist=True)
    elif case == 'case1':
        # line 1: for degenerate divisor Q in J(Fp)
        lcE, const_mult_line, mult_line, sq_line = line_case1(Q_prec, Q, line, c_vec)
    elif case == 'case2' and twist is True:
        # line 2: for general divisor Q in Jt(Fp^s)
        lcE, const_mult_line, mult_line, sq_line = line_case2(Q_prec, Q, line, c_vec, twist=True)
    else:
        # line 2: for general divisor Q in J(Fp)
        lcE, const_mult_line, mult_line, sq_line = line_case2(Q_prec, Q, line, c_vec)

    return R, lcE, const_mult_line, mult_line, sq_line, mult_, sq_


def DBL(T, Q_prec, Q, F, c_vec=None, case: str = 'case1', twist: bool = False):
    """
    :param T: divisor in Fan et al. coordinate system T = [U11, U10, V11, V10, Z11, Z12, z11, z12]
    :param Q_prec: precomputation vector related to divisor Q
    :param Q: degenerate or general divisor for line evaluation
    :param F: curve coefficients C/Fp: y^2 = f(x)
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param twist: twist = True => P in Jt(Fp^s), Q in J(Fp) or twist = False => P in J(Fp), Q in Jt(Fp^s)
    :return: R = [2]*T (in Fan et al. coordinate system) and lcE: line evaluated at Q

    cost: 38 mult. + 6 sq.
    """
    f0, f1, f2, f3, f4 = F[0], F[1], F[2], F[3], F[4]
    U11, U10, V11, V10 = T[0], T[1], T[2], T[3]
    Z11, Z12, z11, z12 = T[4], T[5], T[6], T[7]
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

    R = [U31, U30, V31, V30, Z31, Z32, z31, z32]  # 2*T, result of the doubling

    # cost for doubling + line computation: 38 mult. + 6 sq.
    mult_, sq_ = 38, 6

    if case == 'case1' and twist is True:
        # line 1: for degenerate divisor Q in Jt(Fp^s)
        lcE, const_mult_line, mult_line, sq_line = line_case1(Q_prec, Q, line, c_vec, twist=True)
    elif case == 'case1':
        # line 1: for degenerate divisor Q in J(Fp)
        lcE, const_mult_line, mult_line, sq_line = line_case1(Q_prec, Q, line, c_vec)
    elif case == 'case2' and twist is True:
        # line 2: for general divisor Q in Jt(Fp^s)
        lcE, const_mult_line, mult_line, sq_line = line_case2(Q_prec, Q, line, c_vec, twist=True)
    else:
        # line 2: for general divisor Q in J(Fp)
        lcE, const_mult_line, mult_line, sq_line = line_case2(Q_prec, Q, line, c_vec)

    return R, lcE, const_mult_line, mult_line, sq_line, mult_, sq_


def new_coordinates(D) -> list:
    """
    :param D: a divisor D in general form D = [u, v] = [x^2 + U1*x + U0, V1*x + V0]
    :return: the divisor D_ in Fan et al. coordinate system D_ = [U1, U0, V1, V0, Z1, Z2, z1, z2]
    """
    U1, U0 = D[0][1], D[0][0]
    V1, V0 = D[1][1], D[1][0]
    Z1, Z2 = (1), (1)
    z1, z2 = (Z1 ** 2), (Z2 ** 2)
    D_ = [U1, U0, V1, V0, Z1, Z2, z1, z2]

    return D_


def HEC_random_point(C):
    """
    :param C: the hyperelliptic curve equation over a finite field F => C/F: y^2 = f(x)
    :return: A random point on the curve P = (xP, yP)
    """
    fx = C.hyperelliptic_polynomials()[0]
    while True:
        xP = fx.base_ring().random_element()
        yP2 = fx(xP)
        if yP2.is_square():
            yP = yP2.sqrt()
            return [xP, yP]


def HEC_distinct_random_points(C, g):
    """
    :param C: the hyperelliptic curve equation over a finite field F => C/F: y^2 = f(x)
    :param g: the genus of the curve (g = 2)
    :return: a list of g = 2 (distinct) points res = [P1, P2]
    """
    res = []
    for i in range(1, g + 1):
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
    :param C: the hyperelliptic curve equation over a finite field F => C/F: y^2 = f(x)
    :return: a random divisor D in general form D = [u, v] = [x^2 + U1*x + U0, V1*x + V0]
    """
    fx = C.hyperelliptic_polynomials()[0]
    R = fx.base_ring()['x']
    (x,) = R._first_ngens(1)
    J = C.jacobian()
    points = HEC_distinct_random_points(C, C.genus())
    u = 1
    for point in points:
        u = u * (x - point[0])
    v = fx.parent().lagrange_polynomial(points)

    return J([u, v])
