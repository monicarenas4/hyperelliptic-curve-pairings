from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div, \
    precomputation_degenerate_div
from pairing_types import twisted_ate_cp8, ate_i


def test_twisted_ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller):
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
    :return:
    """
    C, Ct, C8 = curves[0], curves[1], curves[2]
    J, Jt, J8 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp8 = fields[0], fields[1]
    c, c2, c3, c4, c5, c6, c7 = c_vec[0], c_vec[1], c_vec[2], c_vec[3], c_vec[4], c_vec[5], c_vec[6]

    P = JC_random_element(C)
    P = h * P  # Force P to have order r
    P = new_coordinates(P)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            Q = HEC_random_point(Ct)
            xQ, yQ = Q[0], Q[1]
            Q = Ct([xQ, yQ])
            Q_prec = precomputation_degenerate_div(Q)
        else:
            Q = JC_random_element(Ct)
            Q = h_ * Q
            Q_prec = precomputation_general_div(Q)

        pairing_value = twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, case)
        print('pairing value = ', pairing_value)
        pairing_value_naf = twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, case, NAF_rep=True)
        print('pairing value NAF = ', pairing_value_naf)

    return None


def test_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, length_miller, family='k16'):
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
    :return:
    """
    C, Ct, C16 = curves[0], curves[1], curves[2]
    J, Jt, J16 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp2, Fq8 = fields[0], fields[1], fields[2]
    c, c2, c3, c4, c5 = c_vec[0], c_vec[1], c_vec[2], c_vec[3], c_vec[4]
    c6, c7, c8, c9, c10 = c_vec[5], c_vec[6], c_vec[7], c_vec[8], c_vec[9]
    c11, c12, c13, c14, c15 = c_vec[10], c_vec[11], c_vec[12], c_vec[13], c_vec[14]

    Q = JC_random_element(Ct)
    Q = h_ * Q  # Force Q to have order r
    Q = new_coordinates(Q)

    pow = (p ** 16 - 1) // r

    cases = ['case1', 'case2']
    # case1 => Degenerate Divisor
    # case2 =>

    for case in cases:
        if case == 'case1':
            P = HEC_random_point(C)
            P_prec = precomputation_degenerate_div(P)
        else:
            P = JC_random_element(C)
            P = h * P
            P_prec = precomputation_general_div(P)

        pairing_value = ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, case, family=family)
        print('pairing value = ', pairing_value)
        pairing_value_naf = ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, case, NAF_rep=True, family=family)
        print('pairing value NAF = ', pairing_value_naf)

    return None
