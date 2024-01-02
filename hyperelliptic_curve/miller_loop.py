from sage.all import Integer
from jacobian_operations import ADD, DBL
from _utils import NAF


def general_miller_function(P, Q, Q_prec, c_vec, F, length_miller, case: str, twist: str = None, P_neg=None):
    """
    :param P:
    :param Q:
    :param Q_prec:
    :param c_vec:
    :param F:
    :param length_miller:
    :param case:
    :param twist:
    :param P_neg:
    :return:
    """
    if P_neg == None:
        length_miller = Integer(length_miller).digits(2)
    else:
        length_miller = NAF(length_miller)

    T, fc = P, 1

    for i in range(len(length_miller) - 2, -1, -1):
        T, lc = DBL(T, Q_prec, Q, F, c_vec, case=case, twist=twist)
        fc = lc * fc ** 2
        if length_miller[i] == 1:
            T, lc = ADD(P, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc
        elif length_miller[i] == -1:
            T, lc = ADD(P_neg, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc

    return fc
