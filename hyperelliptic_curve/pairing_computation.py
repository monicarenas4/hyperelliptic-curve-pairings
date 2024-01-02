from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div
from pairing_types import Twisted_Ate_cp8, Twisted_Ate_naf_cp8, Ate_i, Ate_i_naf

def test_twisted_ate_cp8(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller):
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
            Q = Ct([xQ,yQ])
            xQ, yQ = xQ / (c2), yQ / (c5)
            Q = C8([xQ,yQ])

            Q = [-xQ, yQ]
            Q_prec = [xQ**2, -xQ**3]
        else:
            Q = JC_random_element(Ct)
            Q = h_ * Q
            u1 = Q[0][1]
            u0 = Q[0][0]
            v1 = Q[1][1]
            v0 = Q[1][0]

            u1 = u1/(c2)
            u0 = u0/(c4)
            v1 = v1/(c3)
            v0 = v0/(c5)
            
            Fp8x = Fp8['x']; (x,) = Fp8x._first_ngens(1)
            ux = x**2 + u1*x + u0
            vx = v1*x + v0

            Q = J8([ux,vx])
            Q_prec = precomputation_general_div(Q)

        pairing_value = Twisted_Ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, case)
        print('pairing value = ', pairing_value)
        pairing_value_naf = Twisted_Ate_naf_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, case)
        print('pairing value NAF = ', pairing_value_naf)
        
    return 0

def test_twisted_ate_k16(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller):
    C, Ct, C16 = curves[0], curves[1], curves[2]
    J, Jt, J16 = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp2, Fq8 = fields[0], fields[1], fields[2]
    c, c2, c3, c4, c5 = c_vec[0], c_vec[1], c_vec[2], c_vec[3], c_vec[4]
    c6, c7, c8, c9, c10 = c_vec[5], c_vec[6], c_vec[7], c_vec[8], c_vec[9]
    c11, c12, c13, c14, c15 = c_vec[10], c_vec[11], c_vec[12], c_vec[13], c_vec[14]

    Q = JC_random_element(Ct)
    Q = h_ * Q  # Force Q to have order r
    Q = new_coordinates(Q)

    pow = (p**16 - 1)//r
    
    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            P = HEC_random_point(C)
            xP, yP = P[0], P[1]
            P = C([xP,yP])
 
            P = [-xP, yP]
            P_prec = [xP**2, -xP**3]
        else:
            P= JC_random_element(C)
            P = h * P
            P_prec = precomputation_general_div(P)

        pairing_value = Ate_i(Q, P, P_prec, c_vec, F, length_miller, U, pow, case)
        print('pairing value = ', pairing_value)
        pairing_value_naf = Ate_i_naf(Q, P, P_prec, c_vec, F, length_miller, U, pow, case)
        print('pairing value NAF = ', pairing_value_naf)

    return 0
