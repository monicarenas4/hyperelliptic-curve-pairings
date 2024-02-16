from write_number_operations import write_number_operations
from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div, \
    precomputation_degenerate_div
from pairing_types import twisted_ate_cp8, ate_i

file_name = 'results/number_of_operations.txt'


def compute_twisted_ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller, k: int):
    """
    :param curves:
    :param jacobians:
    :param fields:
    :param c_vec:
    :param F:
    :param U:
    :param p:
    :param r:
    :param h:
    :param h_:
    :param length_miller:
    :param k: embedding degree
    :return:
    """
    C, Ct, C8 = curves[0], curves[1], curves[2]
    J, Jt, J8 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp8 = fields[0], fields[1]

    P = JC_random_element(C)
    P = h * P  # Force P to have order r
    P = new_coordinates(P)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            Q = HEC_random_point(Ct)
            xQ, yQ = Q[0], Q[1]
            Q_prec, mult_pre, sq_pre = precomputation_degenerate_div(Q)
        else:
            Q = JC_random_element(Ct)
            Q = h_ * Q
            Q_prec, mult_pre, sq_pre = precomputation_general_div(Q)

        write_number_operations(file_name, embedding_degree=k, function='twisted_ate',
                                case=case, mult_pre=mult_pre, sq_pre=sq_pre, total=(mult_pre + sq_pre))

        pairing_value = twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, k=k, case=case)
        pairing_value_naf = twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, k=k, case=case, NAF_rep=True)
        print('pairing value = ', pairing_value)
        print('pairing value NAF = ', pairing_value_naf)

    return None


def compute_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, length_miller, k: int, family='k16'):
    """
    :param curves:
    :param jacobians:
    :param fields:
    :param c_vec:
    :param F:
    :param U:
    :param W:
    :param p:
    :param r:
    :param h:
    :param h_:
    :param length_miller:
    :param k: embedding degree
    :param family: polynomial family
    :return:
    """
    C, Ct, C16 = curves[0], curves[1], curves[2]
    J, Jt, J16 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp2, Fq8 = fields[0], fields[1], fields[2]

    Q = JC_random_element(Ct)
    Q = h_ * Q  # Force Q to have order r
    Q = new_coordinates(Q)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            P = HEC_random_point(C)
            P_prec, mult_pre, sq_pre = precomputation_degenerate_div(P)
        else:
            P = JC_random_element(C)
            P = h * P
            P_prec, mult_pre, sq_pre = precomputation_general_div(P)

        write_number_operations(file_name, embedding_degree=k, function='ate_i',
                                case=case, mult_pre=mult_pre, sq_pre=sq_pre, total=mult_pre + sq_pre)

        pairing_value = ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, k=k, case=case, family=family)
        pairing_value_naf = ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, k=k,
                                  case=case, NAF_rep=True, family=family)
        print('pairing value = ', pairing_value)
        print('pairing value NAF = ', pairing_value_naf)

    return None
