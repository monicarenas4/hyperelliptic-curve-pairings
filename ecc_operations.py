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
    v = xQ - xR

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

    print("They are the same:", xR) if (xR == xP) else None

    lambda_add = (yR - yP) / (xR - xP)
    l = yQ - yR - lambda_add * (xQ - xR)
    R = R + P
    v = xQ - xR

    return R, l, v
