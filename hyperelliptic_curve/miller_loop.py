from sage.all import Integer
from math import log2, floor

from jacobian_operations import ADD, DBL
from _utils import NAF, hamming_weight, NAf_hamming_weight
from write_number_operations import operations_miller_loop, operations_miller_loop_head

file_name = 'results/number_of_operations.txt'


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
    global P_neg, mult_line, sq_line, mult_DBL, sq_DBL, mult_ADD, sq_ADD

    if not NAF_rep:
        vector_miller = Integer(length_miller).digits(2)
    else:
        vector_miller = NAF(length_miller)
        P_neg = [P[0], P[1], -P[2], -P[3], P[4], P[5], P[6], P[7]]

    T, fc = P, 1

    for i in range(len(vector_miller) - 2, -1, -1):
        T, lc, mult_line, sq_line, mult_DBL, sq_DBL = DBL(T, Q_prec, Q, F, c_vec, case=case, twist=twist)
        fc = lc * fc ** 2  # M1 S1
        if vector_miller[i] == 1:
            T, lc, mult_line, sq_line, mult_ADD, sq_ADD = ADD(P, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc  # M2
        elif vector_miller[i] == -1:
            T, lc, mult_line, sq_line, mult_ADD, sq_ADD = ADD(P_neg, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc  # M2

    operations_miller_loop_head(file_name)

    if not NAF_rep:
        operations_miller_loop(file_name, 'miller', case, twist, NAF_rep,
                               floor(log2(length_miller) - 1), hamming_weight(vector_miller),
                               mult_line, sq_line, mult_DBL, sq_DBL, mult_ADD, sq_ADD, 2, 1)
    else:
        operations_miller_loop(file_name, 'miller', case, twist, NAF_rep,
                               floor(log2(length_miller) - 1), NAf_hamming_weight(vector_miller),
                               mult_line, sq_line, mult_DBL, sq_DBL, mult_ADD, sq_ADD, 2, 1)

    return fc
