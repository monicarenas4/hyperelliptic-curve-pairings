def DBL_step(Q, R, a: int = 0):
    """
    :param Q: point in Q
    :param R: point in R
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
    :param P: point in P
    :param Q: point in Q
    :param R: point in R
    :return: l, v
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
    xQ, yQ = Q[0], Q[1]
    xR, yR, zR = R[0], R[1], R[2]
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
    # h = l / v
    return R, l, v


def ADD_miller_point(P, Q, R):
    xP, yP, zP = P[0], P[1], P[2]
    xQ, yQ = Q[0], Q[1]
    xR, yR = R[0], R[1]
    T1 = zP ** 2
    B = xR * T1
    R2 = yR ** 2
    D = ((yR + zP) ** 2 - R2 - T1) * T1
    H = B - xP
    I = H ** 2
    E = 4 * I
    J = H * E
    r = D - 2 * yP
    v = xP * E
    X3 = r ** 2 - J - 2 * v
    Y3 = r * (v - X3) - 2 * yP * J
    Z3 = ((zP + H) ** 2 - T1 - I)
    T3 = Z3 ** 2
    R = [X3, Y3, Z3]
    l = (yQ - yR) * Z3 - (xQ - xR) * r * 1
    v = (xQ * T3 - X3)
    # h = l / v

    return R, l, v
