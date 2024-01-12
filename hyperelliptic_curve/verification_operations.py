from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div, \
    precomputation_degenerate_div
from pairing_types import twisted_ate_cp8, ate_i
from write_number_operations import write_number_operations

from random import randint

file_name = 'results/number_of_operations.txt'


def test_bilinearity_Twisted_Ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller):
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

    P2 = JC_random_element(C)
    P2 = h * P2  # Force P to have order r

    randint_P = randint(0, r - 1)
    P1 = randint_P * P2
    P2 = new_coordinates(P2)
    P1 = new_coordinates(P1)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            Q = HEC_random_point(Ct)
            Q1 = Q
            Q1_prec, mult_pre, sq_pre = precomputation_degenerate_div(Q1)
            Q2 = Q1
            Q2_prec = Q1_prec
            randint_Q = 1
        else:
            Q = JC_random_element(Ct)
            Q2 = h_ * Q
            Q2_prec, mult_pre, sq_pre = precomputation_general_div(Q2)

            randint_Q = randint(0, r - 1)
            Q1 = randint_Q * Q2
            Q1_prec, mult_pre, sq_pre = precomputation_general_div(Q1)

        write_number_operations(file_name, 'test_bilinearity_Twisted_Ate', case, mult_pre=mult_pre, sq_pre=sq_pre)

        print('----------------\nCASE: %s\n----------------' % case)
        pairing_value1 = twisted_ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, case)
        # print('pairing value 1 =\n', pairing_value1)
        pairing_value2 = twisted_ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, case) ** (
                randint_P * randint_Q)
        # print('pairing value 2 =\n', pairing_value2)
        print('bilinearity test:',
              pairing_value1 == pairing_value2) if pairing_value1 != 1 and pairing_value2 != 1 else print('review code')

        pairing_value_naf1 = twisted_ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, case,
                                             NAF_rep=True)
        # print('pairing value NAF 1 =\n', pairing_value_naf1)
        pairing_value_naf2 = twisted_ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, case,
                                             NAF_rep=True) ** (randint_P * randint_Q)
        # print('pairing value NAF 2 =\n', pairing_value_naf2)
        print('bilinearity NAF test:',
              pairing_value_naf1 == pairing_value_naf2) if pairing_value_naf1 != 1 and pairing_value_naf2 != 1 \
            else print('review code')

    return None


def test_bilinearity_Ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, length_miller):
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
    :return:
    """
    C, Ct, C16 = curves[0], curves[1], curves[2]
    J, Jt, J16 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp2, Fq8 = fields[0], fields[1], fields[2]
    c, c2, c3, c4, c5 = c_vec[0], c_vec[1], c_vec[2], c_vec[3], c_vec[4]
    c6, c7, c8, c9, c10 = c_vec[5], c_vec[6], c_vec[7], c_vec[8], c_vec[9]
    c11, c12, c13, c14, c15 = c_vec[10], c_vec[11], c_vec[12], c_vec[13], c_vec[14]

    Q2 = JC_random_element(Ct)
    Q2 = h_ * Q2  # Force P to have order r

    randint_Q = randint(0, r - 1)
    Q1 = randint_Q * Q2
    Q2 = new_coordinates(Q2)
    Q1 = new_coordinates(Q1)

    # pow = (p ** 16 - 1) // r

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            P = HEC_random_point(C)
            xP, yP = P[0], P[1]
            P = C([xP, yP])

            P1 = P
            P1_prec, mult_pre, sq_pre = precomputation_degenerate_div(P1)
            P2 = P1
            P2_prec = P1_prec
            randint_P = 1
        else:
            P = JC_random_element(C)
            P2 = h * P
            P2_prec, mult_pre, sq_pre = precomputation_general_div(P2)

            randint_P = randint(0, r - 1)
            P1 = randint_P * P2
            P1_prec, mult_pre, sq_pre = precomputation_general_div(P1)

        write_number_operations(file_name, 'test_bilinearity_Twisted_Ate', case, mult_pre=mult_pre, sq_pre=sq_pre)


        print('----------------\nCASE: %s\n----------------' % case)
        pairing_value1 = ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, W, case)
        # print('pairing value 1 =\n', pairing_value1)
        pairing_value2 = ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, W, case) ** (
                randint_Q * randint_P)
        # print('pairing value 2 =\n', pairing_value2)
        print('bilinearity test:',
              pairing_value1 == pairing_value2) if pairing_value1 != 1 and pairing_value2 != 1 else print('review code')

        pairing_value_naf1 = ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, W, case, NAF_rep=True)
        # print('pairing value NAF 1 =\n', pairing_value_naf1)
        pairing_value_naf2 = ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, W, case, NAF_rep=True) ** (
                randint_Q * randint_P)
        # print('pairing value NAF 2 =\n', pairing_value_naf2)
        print('bilinearity NAF test:',
              pairing_value_naf1 == pairing_value_naf2) if pairing_value_naf1 != 1 and pairing_value_naf2 != 1 \
            else print('review code')

    return None
