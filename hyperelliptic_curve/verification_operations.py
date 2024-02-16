from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div, \
    precomputation_degenerate_div
from pairing_types import twisted_ate_cp8, ate_i

from random import randint


def test_bilinearity_Twisted_Ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller, k: int):
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
    :param k: embeding degree
    :return:
    """
    C, Ct, C8 = curves[0], curves[1], curves[2]
    J, Jt, J8 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp8 = fields[0], fields[1]

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
            Q2_prec, _, _ = precomputation_general_div(Q2)

            randint_Q = randint(0, r - 1)
            Q1 = randint_Q * Q2
            Q1_prec, _, _ = precomputation_general_div(Q1)

        print('----------------\nCASE: %s\n----------------' % case)
        pairing_value1 = twisted_ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, k=k, case=case)
        pairing_value2 = twisted_ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, k=k,
                                         case=case) ** (randint_P * randint_Q)
        print('bilinearity test:',
              pairing_value1 == pairing_value2) if (pairing_value1 and pairing_value2) != 1 else print('review code')

        pairing_value_naf1 = twisted_ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, k=k,
                                             case=case, NAF_rep=True)
        pairing_value_naf2 = twisted_ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, k=k,
                                             case=case, NAF_rep=True) ** (randint_P * randint_Q)
        print('bilinearity NAF test:',
              pairing_value_naf1 == pairing_value_naf2) if (pairing_value_naf1 and pairing_value_naf2) != 1 \
            else print('review code')

    return None


def test_bilinearity_Ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, length_miller, k: int, family='k16'):
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
    :param k: embeding degree
    :param family:
    :return:
    """
    C, Ct, C16 = curves[0], curves[1], curves[2]
    J, Jt, J16 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp2, Fq8 = fields[0], fields[1], fields[2]

    Q2 = JC_random_element(Ct)
    Q2 = h_ * Q2  # Force P to have order r

    randint_Q = randint(0, r - 1)
    Q1 = randint_Q * Q2
    Q2 = new_coordinates(Q2)
    Q1 = new_coordinates(Q1)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            P = HEC_random_point(C)
            xP, yP = P[0], P[1]
            P = C([xP, yP])
            P1 = P
            P1_prec, _, _ = precomputation_degenerate_div(P1)
            P2 = P1
            P2_prec = P1_prec
            randint_P = 1
        else:
            P = JC_random_element(C)
            P2 = h * P
            P2_prec, _, _ = precomputation_general_div(P2)
            randint_P = randint(0, r - 1)
            P1 = randint_P * P2
            P1_prec, _, _ = precomputation_general_div(P1)

        print('----------------\nCASE: %s\n----------------' % case)

        pairing_value1 = ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, W, k=k, case=case, family=family)
        pairing_value2 = ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, W, k=k, case=case,
                               family=family) ** (randint_Q * randint_P)
        print('bilinearity test:',
              pairing_value1 == pairing_value2) if (pairing_value1 and pairing_value2) != 1 else print('review code')

        pairing_value_naf1 = ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, W, k=k, case=case, NAF_rep=True,
                                   family=family)
        pairing_value_naf2 = ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, W, k=k, case=case, NAF_rep=True,
                                   family=family) ** (randint_Q * randint_P)
        print('bilinearity NAF test:',
              pairing_value_naf1 == pairing_value_naf2) if pairing_value_naf1 != 1 and pairing_value_naf2 != 1 \
            else print('review code')

    return None
