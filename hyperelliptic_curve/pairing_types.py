from final_exponentiation import final_exponentiation_cp8, final_exponentiation_k16
from miller_loop import general_miller_function


def twisted_ate_cp8_general(P, Q, Q_prec, c_vec, F, length_miller, U, K, case: str = 'case1', NAF_rep=False):
    if NAF_rep == False:
        miller_function = general_miller_function(P, Q, Q_prec, c_vec, F, length_miller, case)
    else:
        P_neg = [P[0], P[1], -P[2], -P[3], P[4], P[5], P[6], P[7]]
        miller_function = general_miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, P_neg=P_neg)

    pairing_value = final_exponentiation_cp8(miller_function, U, K)

    return pairing_value


def ate_i_general(Q, P, P_prec, c_vec, F, length_miller, U, K, case: str = 'case1', NAF_rep=False):
    if NAF_rep == False:
        miller_function = general_miller_function(Q, P, P_prec, c_vec, F, length_miller, case, twist='k16')
    else:
        Q_neg = [Q[0], Q[1], -Q[2], -Q[3], Q[4], Q[5], Q[6], Q[7]]
        miller_function = general_miller_function(Q, P, P_prec, c_vec, F, length_miller, case, twist='k16', P_neg=Q_neg)

    pairing_value = final_exponentiation_k16(miller_function, U, K)

    return pairing_value
