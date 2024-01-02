from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div
from pairing_types import Twisted_Ate_cp8, Twisted_Ate_naf_cp8, Ate_i, Ate_i_naf
from random import randint


def test_bilinearity_Twisted_Ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller):
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
            xQ, yQ = Q[0], Q[1]
            Q = Ct([xQ, yQ])
            xQ, yQ = xQ / (c2), yQ / (c5)
            Q = C8([xQ, yQ])

            Q1 = [-xQ, yQ]
            Q1_prec = [xQ ** 2, -xQ ** 3]
            Q2 = Q1
            Q2_prec = Q1_prec
            randint_Q = 1
        else:
            Q = JC_random_element(Ct)
            u1 = Q[0][1]
            u0 = Q[0][0]
            v1 = Q[1][1]
            v0 = Q[1][0]

            u1 = u1 / (c2)
            u0 = u0 / (c4)
            v1 = v1 / (c3)
            v0 = v0 / (c5)

            Fp8x = Fp8['x']
            (x,) = Fp8x._first_ngens(1)
            ux = x ** 2 + u1 * x + u0
            vx = v1 * x + v0

            Q2 = J8([ux, vx])
            Q2 = h_ * Q2
            Q2_prec = precomputation_general_div(Q2)

            randint_Q = randint(0, r - 1)
            Q1 = randint_Q * Q2
            Q1_prec = precomputation_general_div(Q1)

        pairing_value1 = Twisted_Ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, case)
        print('pairing value 1 = ', pairing_value1)
        pairing_value2 = Twisted_Ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, case) ** (
                    randint_P * randint_Q)
        print('pairing value 2 = ', pairing_value2)
        print('bilinearity test: ', pairing_value1 == pairing_value2)

        pairing_value_naf1 = Twisted_Ate_naf_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, case)
        print('pairing value NAF 1 = ', pairing_value_naf1)
        pairing_value_naf2 = Twisted_Ate_naf_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, case) ** (
                    randint_P * randint_Q)
        print('pairing value NAF 2 = ', pairing_value_naf2)
        print('bilinearity NAF test: ', pairing_value_naf1 == pairing_value_naf2)

    return 0


def test_bilinearity_Ate_i(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller):
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

    pow = (p ** 16 - 1) // r

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            P = HEC_random_point(C)
            xP, yP = P[0], P[1]
            P = C([xP, yP])

            P1 = [-xP, yP]
            P1_prec = [xP ** 2, -xP ** 3]
            P2 = P1
            P2_prec = P1_prec
            randint_P = 1
        else:
            P = JC_random_element(C)
            P2 = h * P
            P2_prec = precomputation_general_div(P2)

            randint_P = randint(0, r - 1)
            P1 = randint_P * P2
            P1_prec = precomputation_general_div(P1)

        pairing_value1 = Ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, pow, case)
        print('pairing value 1 = ', pairing_value1)
        pairing_value2 = Ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, pow, case) ** (randint_Q * randint_P)
        print('pairing value 2 = ', pairing_value2)
        print('bilinearity test: ', pairing_value1 == pairing_value2)

        pairing_value_naf1 = Ate_i_naf(Q1, P1, P1_prec, c_vec, F, length_miller, U, pow, case)
        print('pairing value NAF 1 = ', pairing_value_naf1)
        pairing_value_naf2 = Ate_i_naf(Q2, P2, P2_prec, c_vec, F, length_miller, U, pow, case) ** (
                    randint_Q * randint_P)
        print('pairing value NAF 2 = ', pairing_value_naf2)
        print('bilinearity NAF test: ', pairing_value_naf1 == pairing_value_naf2)

    return 0
