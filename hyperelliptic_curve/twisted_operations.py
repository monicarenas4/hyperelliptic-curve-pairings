from math import floor, log2
from operations import ADD, DBL
from final_exponentiation import final_exponentiation_k8
from final_exponentiation import final_exponentiation_k16
from final_exponentiation import final_exponentiation_k16_FK


def twisted_ate_k8(D, Dn, Q_vec, Q, F, s, U, K, L, case: str = 'case1'):
    """
    :param D:
    :param Dn:
    :param Q_vec:
    :param Q:
    :param F:
    :param s:
    :param U:
    :param K:
    :param L:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :return:
    """
    T, fc = D, 1

    for i in range(len(s) - 2, -1, -1):
        T, lc = DBL(T, Q_vec, Q, F, L, case)
        fc = lc * fc ** 2
        if s[i] == 1:
            T, lc = ADD(D, T, Q_vec, Q, F, L, case)
            fc = lc * fc
        elif s[i] == -1:
            T, lc = ADD(Dn, T, Q_vec, Q, F, L, case)
            fc = lc * fc

    fc = final_exponentiation_k8(fc, U, K)

    return fc


def twisted_ate_k16(D, E, F, ML, L, U, K):
    """
    :param D:
    :param E:
    :param F:
    :param ML:
    :param L:
    :param U:
    :param K:
    :return:
    """
    T, fc = D, 1
    s = [int(x) for x in "{0:0b}".format(ML)]

    for i in range(len(s) - 2, -1, -1):
        T, lc = DBL(T, E, L, F, L)
        fc = lc * fc ** 2
        if s[i] == 1:
            T, lc = ADD(D, T, E, L, F, L)
            fc = lc * fc

    fc = final_exponentiation_k16(fc, U, K)

    return fc


def twisted_ate_k16_FK(D, Dn, E, F, ML, L, s, U, K):
    """
    :param D:
    :param Dn:
    :param E:
    :param F:
    :param ML:
    :param L:
    :param s:
    :param U:
    :param K:
    :return:
    """
    T, fc = D, 1

    for i in range(len(s) - 2, -1, -1):
        T, lc = DBL(T, E, L, F)
        fc = lc * fc ** 2
        if s[i] == 1:
            T, lc = ADD(D, T, E, L, F)
            fc = lc * fc
        elif s[i] == -1:
            T, lc = ADD(Dn, T, E, L, F)
            fc = lc * fc

    fc = final_exponentiation_k16_FK(fc, U, K)

    return fc
