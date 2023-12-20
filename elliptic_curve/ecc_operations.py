def DBL_step(Q, R, a: int = 0):
    """
    :param Q: Q-point with (x,y,z)-coordinates
    :param R: R-point with (x,y,z)-coordinates
    :param a: intercept
    :return: R, l, v
    """
    xQ, yQ = Q[0], Q[1]
    xR, yR = R[0], R[1]

    lambda_dbl = (3 * xR ** 2 + a) / (2 * yR)
    l = yQ - yR - lambda_dbl * (xQ - xR)
    R = 2 * R
    v = xQ - R[0]

    return R, l, v


def ADD_step(P, Q, R):
    """
    :param P: P-point with (x,y,z)-coordinates
    :param Q: Q-point with (x,y,z)-coordinates
    :param R: R-point with (x,y,z)-coordinates
    :return: R, l, v
    """
    xP, yP = P[0], P[1]
    xR, yR = R[0], R[1]
    xQ, yQ = Q[0], Q[1]

    if xR != xP:
        lambda_add = (yR - yP) / (xR - xP)
        l = yQ - yR - lambda_add * (xQ - xR)
        R = R + P
        v = xQ - R[0]
    else:
        l = xQ - xP
        R = R + P
        v = 1

    return R, l, v


def DBL_miller_point(Q, R, a: int = 0):
    """
    :param Q: Q-point with (x,y,z)-coordinates
    :param R: R-point with (x,y,z)-coordinates
    :param a: int
    :return: R, l, v
    """
    xQ, yQ, _ = Q
    xR, yR, zR = R
    T1 = zR ** 2
    A = yR ** 2
    B = xR * A
    C = 3 * xR ** 2 + a * T1 ** 2
    X3 = C ** 2 - 8 * B
    Z3 = (yR + zR) ** 2 - A - T1
    Y3 = C * (4 * B - X3) - 8 * A ** 2
    T3 = Z3 ** 2
    R = [X3, Y3, Z3]
    l = T3 * (yQ * Z3 - C * xQ) + Y3 + C * X3
    v = (xQ * T3 - X3)

    return R, l, v


def ADD_miller_point(P, Q, R):
    """
    :param P: P-point with (x,y,z)-coordinates
    :param Q: Q-point with (x,y,z)-coordinates
    :param R: R-point with (x,y,z)-coordinates
    :return: R, l, v
    """
    xP, yP, _ = P
    xQ, yQ, _ = Q
    xR, yR, zR = R
    T1 = zR ** 2
    B = xP * T1
    R2 = yP ** 2
    D = ((yP + zR) ** 2 - R2 - T1) * T1
    H = B - xR
    I = H ** 2
    E = 4 * I
    J = H * E
    r = D - 2 * yR
    v = xR * E
    X3 = r ** 2 - J - 2 * v
    Y3 = r * (v - X3) - 2 * yR * J
    Z3 = ((zR + H) ** 2 - T1 - I)
    T3 = Z3 ** 2
    R = [X3, Y3, Z3]
    l = (yQ - yP) * Z3 - (xQ - xP) * r * 1
    v = (xQ * T3 - X3)

    return R, l, v
