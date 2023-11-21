def DBL_step(R, Q, a: int):
    """
    R: point in R
    Q: point in Q
    a: intercept
    """
    xR, yR = R[0], R[1]
    xQ, yQ = Q[0], Q[1]

    lambda_dbl = (3 * xR + a) / (2 * yR)
    l = yQ - yR - lambda_dbl * (xQ - xR)
    R = 2 * R
    v = xQ - xR

    return R, l, v


def ADD_step(P, R, Q):
    """
    P: point in P
    R: point in R
    Q: point in Q
    """
    xP, yP = P[0], P[1]
    xR, yR = R[0], R[1]
    xQ, yQ = Q[0], Q[1]

    lambda_add = (yR - yP) / (xR - xP)
    l = yQ - yR - lambda_add * (xQ - xR)
    R = R + P
    v = xQ - R[0]

    return l, v
