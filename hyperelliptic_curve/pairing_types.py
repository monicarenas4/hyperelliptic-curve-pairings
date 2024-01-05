from final_exponentiation import final_exponentiation_cp8, final_exponentiation_k16
from miller_loop import miller_function


def twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, K, case: str = 'case1', NAF_rep=False):
    """
    :param P:
    :param Q:
    :param Q_prec:
    :param c_vec:
    :param F:
    :param length_miller:
    :param U:
    :param K:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param NAF_rep: true or false
    :return:
    """
    if not NAF_rep:
        miller_fun = miller_function(P, Q, Q_prec, c_vec, F, length_miller, case)
    else:
        miller_fun = miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, NAF_rep=NAF_rep)

    pairing_value = final_exponentiation_cp8(miller_fun, U, K)

    return pairing_value


def ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, case: str = 'case1', NAF_rep=False):
    """
    :param Q:
    :param P:
    :param P_prec:
    :param c_vec:
    :param F:
    :param length_miller:
    :param U:
    :param K:
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param NAF_rep: true or false
    :return:
    """
    if NAF_rep == False:
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case, twist='k16')
    else:
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case, twist='k16', NAF_rep=NAF_rep)

    pairing_value = final_exponentiation_k16(miller_fun, U, W)

    return pairing_value
