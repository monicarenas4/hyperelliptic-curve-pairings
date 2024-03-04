from final_exponentiation import final_exponentiation_cp8
from final_exponentiation import final_exponentiation_k16, final_exponentiation_new_k16
from final_exponentiation import final_exponentiation_k24

from miller_loop import miller_function
from write_number_operations import write_number_operations

file_name = 'results/number_of_operations.txt'


def twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, k: int, case: str = 'case1', NAF_rep=False):
    """
    :param P:
    :param Q:
    :param Q_prec:
    :param c_vec: powers of C
    :param F: coefficients of C
    :param length_miller:
    :param U: related to the seed u
    :param Fp: field
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param NAF_rep: true or false
    :return:
    """
    if not NAF_rep:
        miller_fun = miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, k)
        pairing_value, exp_u, exp_u0, exp_lx, mult, sq, inv, frob, total = final_exponentiation_cp8(miller_fun, U, Fp)
    else:
        miller_fun = miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, k, NAF_rep=NAF_rep)
        pairing_value, exp_u, exp_u0, exp_lx, mult, sq, inv, frob, total = final_exponentiation_cp8(miller_fun, U, Fp,
                                                                                                    NAF_rep=NAF_rep)

    write_number_operations(file_name, function='final_exp twisted_ate_cp8', embedding_degree=k, case=case,
                            NAF_rep=NAF_rep, exp_u=exp_u, exp_u0=exp_u0, exp_lx=exp_lx, mult_FE=mult, sq_FE=sq,
                            inv_FE=inv, frobenius=frob, total=total)

    return pairing_value


def ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, k: int, case: str = 'case1', NAF_rep=False, family='k16'):
    """
    :param Q:
    :param P:
    :param P_prec:
    :param c_vec:
    :param F:
    :param length_miller:
    :param U:
    :param W:
    :param k: embedding degree
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param NAF_rep: true or false
    :param family:
    :return:
    """
    global pairing_value

    if family == "k16" and NAF_rep == True:
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True,
                                     NAF_rep=NAF_rep)
        pairing_value, exp_u, exp_um, mult, sq, inv, frob_power, total = final_exponentiation_k16(miller_fun, U, W,
                                                                                                  NAF_rep=NAF_rep)
        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, NAF_rep=NAF_rep,
                                exp_u=exp_u, exp_um=exp_um, mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power,
                                total=total)

    elif family == "k16":
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True)
        pairing_value, exp_u, exp_um, mult, sq, inv, frob_power, total = final_exponentiation_k16(miller_fun, U, W)
        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u, exp_um=exp_um,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "new_k16" and NAF_rep == True:
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case,
                                     k=k, twist=True, NAF_rep=NAF_rep)
        pairing_value, exp_u, exp_up, mult, sq, inv, frob_power, total = final_exponentiation_new_k16(miller_fun, U, W,
                                                                                                      NAF_rep=NAF_rep)
        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case,
                                NAF_rep=NAF_rep, exp_u=exp_u, exp_up=exp_up, mult_FE=mult,
                                sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "new_k16":
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True)
        pairing_value, exp_u, exp_up, mult, sq, inv, frob_power, total = final_exponentiation_new_k16(miller_fun, U, W)

        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u, exp_up=exp_up,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "k24" and NAF_rep == True:
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True, NAF_rep=NAF_rep)
        pairing_value, exp_u, mult, sq, inv, frob_power, total = final_exponentiation_k24(miller_fun, U, W)

        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "k24":
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True)
        pairing_value, exp_u, mult, sq, inv, frob_power, total = final_exponentiation_k24(miller_fun, U, W)

        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    return pairing_value
