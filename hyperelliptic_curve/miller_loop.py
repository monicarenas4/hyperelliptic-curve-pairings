from sage.all import Integer
from math import log2, floor

from jacobian_operations import ADD, DBL
from _utils import NAF, hamming_weight, NAf_hamming_weight
from write_number_operations import write_number_operations

file_name = 'results/number_of_operations.txt'


def miller_function(P, Q, Q_prec, c_vec, F, length_miller, case: str, k: int, twist: str = None, NAF_rep=False):
    """
    :param P: point over Fp
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
    global P_neg, mult_DBL, sq_DBL, mult_ADD, sq_ADD, mult_line_DBL, sq_line_DBL
    global mult_line_ADD, sq_line_ADD, sk, mk, mult_miller_DBL

    if not NAF_rep:
        vector_miller = Integer(length_miller).digits(2)
    else:
        vector_miller = NAF(length_miller)
        P_neg = [P[0], P[1], -P[2], -P[3], P[4], P[5], P[6], P[7]]

    T, fc = P, 1
    number_of_DBL = floor(log2(length_miller)) - 1

    for i in range(len(vector_miller) - 2, -1, -1):
        T, lc, mult_line_DBL, sq_line_DBL, mult_DBL, sq_DBL = DBL(T, Q_prec, Q, F, c_vec, case=case, twist=twist)
        fc = lc * fc ** 2  # 1m_8, 1s_8 (sparse)
        if vector_miller[i] == 1:
            T, lc, mult_line_ADD, sq_line_ADD, mult_ADD, sq_ADD = ADD(P, T, Q_prec, Q, F, c_vec, case=case, twist=twist)
            fc = lc * fc  # 1m_8 (sparse)
        elif vector_miller[i] == -1:
            T, lc, mult_line_ADD, sq_line_ADD, mult_ADD, sq_ADD = ADD(P_neg, T, Q_prec, Q, F, c_vec, case=case,
                                                                      twist=twist)
            fc = lc * fc  # 1m_8

    if not NAF_rep:
        number_of_ADD = hamming_weight(vector_miller) - 1
    else:
        number_of_ADD = NAf_hamming_weight(vector_miller) - 1

    if k == 8:
        mk, mk_DBL, sk, fk, sk_cyclo, ik = 27, 18, 18, 6, 12, 69  # check
        mult_miller_DBL = (mk_DBL * (number_of_DBL - 1))  # check
    elif k == 16:
        mk, sk, fk, sk_cyclo, ik = 81, 54, 14, 36, 159
        mult_miller_DBL = (mk * (number_of_DBL - 1))

    mult_DBL = mult_DBL * number_of_DBL
    sq_DBL = sq_DBL * number_of_DBL
    mult_ADD = mult_ADD * number_of_ADD
    sq_ADD = sq_ADD * number_of_ADD
    # mult_miller_DBL = (mk * (number_of_DBL - 1))  # check
    sq_miller_DBL = (sk * (number_of_DBL - 1))
    mult_miller_ADD = (mk * number_of_ADD)
    sq_miller_ADD = 0 * number_of_ADD

    total = mult_DBL + sq_DBL + mult_ADD + sq_ADD + mult_miller_DBL + sq_miller_DBL + mult_miller_ADD + sq_miller_ADD

    write_number_operations(file_name,
                            embedding_degree=k,
                            function='miller', case=case, twist=twist, NAF_rep=NAF_rep,
                            number_of_DBL=number_of_DBL,
                            number_of_ADD=number_of_ADD,
                            mult_DBL=mult_DBL, sq_DBL=sq_DBL,
                            mult_ADD=mult_ADD, sq_ADD=sq_ADD,
                            mult_miller_DBL=mult_miller_DBL,
                            sq_miller_DBL=sq_miller_DBL,
                            mult_miller_ADD=mult_miller_ADD,
                            sq_miller_ADD=sq_miller_ADD, total=total)

    return fc
