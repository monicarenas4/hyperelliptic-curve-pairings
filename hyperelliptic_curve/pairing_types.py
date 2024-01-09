from final_exponentiation import final_exponentiation_cp8, final_exponentiation_k16
from miller_loop import miller_function
from write_number_operations import operations_bilinearity_check, operations_bilinearity_check_head

file_name = 'results/number_of_operations.txt'


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

    pairing_value, mult_FE, sq_FE = final_exponentiation_cp8(miller_fun, U, K)

    operations_bilinearity_check_head(file_name)
    operations_bilinearity_check(file_name, 'final_exp twisted_ate_cp8', case, NAF_rep, mult_FE, sq_FE)

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

    pairing_value, mult_FE, sq_FE = final_exponentiation_k16(miller_fun, U, W)

    operations_bilinearity_check_head(file_name)
    operations_bilinearity_check(file_name, 'final_exp ate_i', case, NAF_rep, mult_FE, sq_FE)

    return pairing_value
