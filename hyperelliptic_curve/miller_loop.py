from sage.all import Integer
from jacobian_operations import ADD, DBL
from _utils import NAF


def miller_function(P, Q, Q_prec, c_vec, F, length_miller, case: str, twist: str = None, NAF_rep=False):
    """
    :param P:
    :param Q:
    :param Q_prec:
    :param c_vec:
    :param F:
    :param length_miller:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param twist:
    :param NAF_rep:
    :return:
    """
    if not NAF_rep:
        length_miller = Integer(length_miller).digits(2)
    else:
        length_miller = NAF(length_miller)

    T, fc = P, 1
    P_neg = [P[0], P[1], -P[2], -P[3], P[4], P[5], P[6], P[7]]

    for i in range(len(length_miller) - 2, -1, -1):
        T, lc = DBL(T, Q_prec, Q, F, c_vec, case=case, twist=twist)
        fc = lc * fc ** 2  # M1 S1
        if length_miller[i] == 1:
            T, lc = ADD(P, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc  # M2
        elif length_miller[i] == -1:
            T, lc = ADD(P_neg, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc  # M2

    return fc
