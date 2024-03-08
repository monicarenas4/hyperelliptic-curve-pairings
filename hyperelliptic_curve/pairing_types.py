from final_exponentiation import final_exponentiation_cp8
from final_exponentiation import final_exponentiation_KT16, final_exponentiation_New16
from final_exponentiation import final_exponentiation_New24

from miller_loop import miller_function
from write_number_operations import write_number_operations

file_name = 'results/number_of_operations.txt'


def twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, k: int, case: str = 'case1', NAF_rep=False):
    """
    :param P: divisor in J(Fp) in Fan et al. coordinate system representation
    :param Q: divisor in Jt(Fp^s) represented as degenerate of general divisor
    :param Q_prec: precomputation values required for line evaluation in Miller loop
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param F: coefficients of the polynomial f(x), s.t. C/Fp: y^2 = f(x)
    :param length_miller: integer that defined the length of the Miller loop. In all cases this is the seed u
    :param U: integers required in the final exponentiation, including the seed u
    :param Fp: base field
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param NAF_rep: true or false
    :return: pairing_value, the output of the pairing computation: e(P, Q) = f_{u,Q}(P)^((p^k - 1)/r)
    """
    if not NAF_rep:
        # compute Miller function f for binary representation of length_miller
        miller_fun = miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, k)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, exp_u0, exp_lx, mult, sq, inv, frob, total = final_exponentiation_cp8(miller_fun, U, Fp)
    else:
        # compute Miller function f for NAF representation of length_miller
        miller_fun = miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, k, NAF_rep=NAF_rep)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, exp_u0, exp_lx, mult, sq, inv, frob, total = final_exponentiation_cp8(miller_fun, U, Fp,
                                                                                                    NAF_rep=NAF_rep)

    write_number_operations(file_name, function='final_exp twisted_ate_cp8', embedding_degree=k, case=case,
                            NAF_rep=NAF_rep, exp_u=exp_u, exp_u0=exp_u0, exp_lx=exp_lx, mult_FE=mult, sq_FE=sq,
                            inv_FE=inv, frobenius=frob, total=total)

    return pairing_value


def ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, k: int, case: str = 'case1', NAF_rep=False, family='KT16'):
    """
    :param Q: divisor in Jt(Fp^s) in Fan et al. coordinate system representation
    :param P: divisor in J(Fp) represented as degenerate of general divisor
    :param P_prec: precomputation values required for line evaluation in Miller loop
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param F: coefficients of the polynomial f(x), s.t. C/Fp: y^2 = f(x)
    :param length_miller: integer that defined the length of the Miller loop. In all cases this is the seed u
    :param U: integers required in the final exponentiation, including the seed u
    :param W: powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    :param k: embedding degree
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param NAF_rep: true or false
    :param family: polynomial family; used to distinguish between different constructions
    :return: pairing_value, the output of the pairing computation: e(Q, P) = f_{u,P}(Q)^((p^k - 1)/r)
    """
    global pairing_value

    if family == "KT16" and NAF_rep == True:
        # compute Miller function f for KT16 in NAF representation of length_miller
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True,
                                     NAF_rep=NAF_rep)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, exp_um, mult, sq, inv, frob_power, total = final_exponentiation_KT16(miller_fun, U, W,
                                                                                                  NAF_rep=NAF_rep)
        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, NAF_rep=NAF_rep,
                                exp_u=exp_u, exp_um=exp_um, mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power,
                                total=total)

    elif family == "KT16":
        # compute Miller function f for KT16 in binary representation of length_miller
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, exp_um, mult, sq, inv, frob_power, total = final_exponentiation_KT16(miller_fun, U, W)
        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u, exp_um=exp_um,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "New16" and NAF_rep == True:
        # compute Miller function f for New16 in NAF representation of length_miller
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case,
                                     k=k, twist=True, NAF_rep=NAF_rep)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, exp_up, mult, sq, inv, frob_power, total = final_exponentiation_New16(miller_fun, U, W,
                                                                                                      NAF_rep=NAF_rep)
        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case,
                                NAF_rep=NAF_rep, exp_u=exp_u, exp_up=exp_up, mult_FE=mult,
                                sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "New16":
        # compute Miller function f for New16 in binary representation of length_miller
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, exp_up, mult, sq, inv, frob_power, total = final_exponentiation_New16(miller_fun, U, W)

        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u, exp_up=exp_up,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "New24" and NAF_rep == True:
        # compute Miller function f for New24 in NAF representation of length_miller
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True, NAF_rep=NAF_rep)
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, mult, sq, inv, frob_power, total = final_exponentiation_New24(miller_fun, U, W)

        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    elif family == "New24":
        # compute Miller function f for New24 in binary representation of length_miller
        miller_fun = miller_function(Q, P, P_prec, c_vec, F, length_miller, case=case, k=k, twist=True) \
        # raise Miller function f to exponent (p^k - 1)/r
        pairing_value, exp_u, mult, sq, inv, frob_power, total = final_exponentiation_New24(miller_fun, U, W)

        write_number_operations(file_name, 'final_exp ate_i', embedding_degree=k, case=case, exp_u=exp_u,
                                mult_FE=mult, sq_FE=sq, inv_FE=inv, frobenius=frob_power, total=total)

    return pairing_value
