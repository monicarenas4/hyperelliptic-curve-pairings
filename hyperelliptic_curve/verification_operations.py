from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div, \
    precomputation_degenerate_div
from pairing_types import twisted_ate_cp8, ate_i

from random import randint

def test_bilinearity_twisted_ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller, k: int):
    """
    :param curves: three curves C/Fp, Ct/Fp^s, Ck/Fp^k, where s = k / d and d = 8 (degree of twist)
    :param jacobians: three Jacobians J/Fp, Jt/Fp^s, Jk/Fp^k, where s = k / d and d = 8 (degree of twist)
    :param fields: three finite field Fp, Fq = Fp^s, Fq^d = Fp^k, where s = k / d and d = 8 (degree of twist)
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param F: coefficients of the polynomial f(x), s.t. C/Fp: y^2 = f(x)
    :param U: integers required in the final exponentiation, including the seed u
    :param p: base field prime
    :param r: prime that divides the order of the Jacobian J/Fp i.e., #J(Fp) = h*r
    :param h: cofactor of the order #J(Fp)
    :param h_: cofactor of the order #J(Fp^k)
    :param length_miller: integer that defined the length of the Miller loop. In all cases this is the seed u
    :param k: embeding degree
    :return:

    verify bilinearity: e([a]*P, [b]*Q) = e(P, Q)^ab
    """
    C, Ct, Ck = curves[0], curves[1], curves[2]
    J, Jt, Jk = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fpk = fields[0], fields[1]

    # Choose random divisor P2 in Jacobian J(Fp)
    P2 = JC_random_element(C)
    P2 = h * P2  # Force P2 to have order r

    randint_P = randint(0, r - 1) # generate random a in Z_r
    P1 = randint_P * P2              # set P1 = [a]*P2
    P2 = new_coordinates(P2)         # Convert P2 to Fan et al. coordinate system P2 = [U21, U20, V21, V20, Z21, Z22, z21, z22]
    P1 = new_coordinates(P1)         # Convert P1 to Fan et al. coordinate system P1 = [U11, U10, V11, V10, Z11, Z12, z11, z12]

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            # second pairing inputs are the same Q1 = Q2
            Q = HEC_random_point(Ct)                                      # random point on twist curve Ct/Fp^s
            Q1 = Q
            Q1_prec, mult_pre, sq_pre = precomputation_degenerate_div(Q1) # precomputation: values required in line evaluation
            Q2 = Q1
            Q2_prec = Q1_prec
            randint_Q = 1
        else:
            # case 2 => General Divisor
            # second pairing inputs are Q2 and Q1 = [b]*Q2
            Q = JC_random_element(Ct)                      # random divisor on twist Jacobian Jt/Fp^s => Q = [x^2 + u1*x + u0, v1*x + v0]
            Q2 = h_ * Q                                    # Force Q2 to have order r
            Q2_prec, _, _ = precomputation_general_div(Q2) # precomputation: values related to Q2 required in line evaluation

            randint_Q = randint(0, r - 1)               # generate random b in Z_r
            Q1 = randint_Q * Q2                            # Q1 = [b]*Q2
            Q1_prec, _, _ = precomputation_general_div(Q1) # precomputation: values related to Q1 required in line evaluation

        print('----------------\nCASE: %s\n----------------' % case)
        # case 1 / case 2: binary representation of length_miller

        # compute e(P1, Q1) = e([a]*P,Q) (case 1) or e(P1, Q1) = e([a]*P,[b]*Q) (case 2)
        pairing_value1 = twisted_ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, k=k, case=case)

        # compute e(P2, Q2) = e(P,Q)^a (case 1) or e(P2, Q2) = e(P,Q)^ab (case 2)
        pairing_value2 = twisted_ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, k=k,
                                         case=case) ** (randint_P * randint_Q)

        # verify: e(P1, Q1) = e(P2, Q2) <=> e([a]*P,Q) = e(P,Q)^a (case1)
        # verify: e(P1, Q1) = e(P2, Q2) <=> e([a]*P,[b]*Q) = e(P,Q)^ab (case2)
        print('bilinearity test:',
              pairing_value1 == pairing_value2) if (pairing_value1 and pairing_value2) != 1 else print('review code')

        # case 1 / case 2: NAF representation of length_miller
        # compute e(P1, Q1) = e([a]*P,Q) (case 1) or e(P1, Q1) = e([a]*P,[b]*Q) (case 2)
        pairing_value_naf1 = twisted_ate_cp8(P1, Q1, Q1_prec, c_vec, F, length_miller, U, Fp, k=k,
                                             case=case, NAF_rep=True)
        # compute e(P2, Q2) = e(P,Q)^a (case 1) or e(P2, Q2) = e(P,Q)^ab (case 2)
        pairing_value_naf2 = twisted_ate_cp8(P2, Q2, Q2_prec, c_vec, F, length_miller, U, Fp, k=k,
                                             case=case, NAF_rep=True) ** (randint_P * randint_Q)
        # verify: e(P1, Q1) = e(P2, Q2) <=> e([a]*P,Q) = e(P,Q)^a (case1)
        # verify: e(P1, Q1) = e(P2, Q2) <=> e([a]*P,[b]*Q) = e(P,Q)^ab (case2)
        print('bilinearity NAF test:',
              pairing_value_naf1 == pairing_value_naf2) if (pairing_value_naf1 and pairing_value_naf2) != 1 \
            else print('review code')

    return None


def test_bilinearity_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, length_miller, k: int, family='KT16'):
    """
    :param curves: three curves C/Fp, Ct/Fp^s, Ck/Fp^k, where s = k / d and d = 8 (degree of twist)
    :param jacobians: three Jacobians J/Fp, Jt/Fp^s, Jk/Fp^k, where s = k / d and d = 8 (degree of twist)
    :param fields: three finite field Fp, Fq = Fp^s, Fq^d = Fp^k, where s = k / d and d = 8 (degree of twist)
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param F: coefficients of the polynomial f(x), s.t. C/Fp: y^2 = f(x)
    :param U: integers required in the final exponentiation, including the seed u
    :param W: powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    :param p: base field prime
    :param r: prime that divides the order of the Jacobian J/Fp i.e., #J(Fp) = h*r
    :param h: cofactor of the order #J(Fp)
    :param h_: cofactor of the order #J(Fp^k)
    :param length_miller: integer that defined the length of the Miller loop. In all cases this is the seed u
    :param k: embeding degree
    :param family: polynomial family; used to distinguish between different constructions
    :return:

    verify bilinearity: e([a]*Q, [b]*P) = e(Q, P)^ab
    """
    C, Ct, Ck = curves[0], curves[1], curves[2]
    J, Jt, Jk = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fp2, Fqk = fields[0], fields[1], fields[2]

    # Choose random divisor Q2 in Jacobian Jt(Fp^s)
    Q2 = JC_random_element(Ct)
    Q2 = h_ * Q2  # Force Q2 to have order r

    randint_Q = randint(0, r - 1) # generate random a in Z_r
    Q1 = randint_Q * Q2              # set Q1 = [a]*Q2
    Q2 = new_coordinates(Q2)         # Convert Q2 to Fan et al. coordinate system Q2 = [U21, U20, V21, V20, Z21, Z22, z21, z22]
    Q1 = new_coordinates(Q1)         # Convert Q1 to Fan et al. coordinate system Q1 = [U11, U10, V11, V10, Z11, Z12, z11, z12]

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 => Degenerate Divisor
            # second pairing inputs are the same P1 = P2
            P = HEC_random_point(C)                           # random point on curve C/Fp
            P1 = P
            P1_prec, _, _ = precomputation_degenerate_div(P1) # precomputation: values related to P1 required in line evaluation
            P2 = P1
            P2_prec = P1_prec
            randint_P = 1
        else:
            # case 2 => General Divisor
            # second pairing inputs are P2 and P1 = [b]*P2
            P = JC_random_element(C)                       # random divisor on Jacobian J/Fp => P = [x^2 + u1*x + u0, v1*x + v0]
            P2 = h * P                                     # Force Q2 to have order r
            P2_prec, _, _ = precomputation_general_div(P2) # precomputation: values related to P2 required in line evaluation
            randint_P = randint(0, r - 1)               # generate random b in Z_r
            P1 = randint_P * P2                            # P1 = [b]*P2
            P1_prec, _, _ = precomputation_general_div(P1) # precomputation: values related to P1 required in line evaluation

        print('----------------\nCASE: %s\n----------------' % case)
        # case 1 / case 2: binary representation of length_miller

        # compute e(Q1, P1) = e([a]*Q,P) (case 1) or e(Q1, P1) = e([a]*Q,[b]*P) (case 2)
        pairing_value1 = ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, W, k=k, case=case, family=family)

        # compute e(Q2, P2) = e(Q,P)^a (case 1) or e(Q2, P2) = e(Q,P)^ab (case 2)
        pairing_value2 = ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, W, k=k, case=case,
                               family=family) ** (randint_Q * randint_P)

        # verify: e(Q1, P1) = e(Q2, P2) <=> e([a]*Q,P) = e(Q,P)^a (case1)
        # verify: e(Q1, P1) = e(Q2, P2) <=> e([a]*Q,[b]*P) = e(Q,P)^ab (case2)
        print('bilinearity test:',
              pairing_value1 == pairing_value2) if (pairing_value1 and pairing_value2) != 1 else print('review code')

        # case 1 / case 2: NAF representation of length_miller
        # compute e(Q1, P1) = e([a]*Q,P) (case 1) or e(Q1, P1) = e([a]*Q,[b]*P) (case 2)
        pairing_value_naf1 = ate_i(Q1, P1, P1_prec, c_vec, F, length_miller, U, W, k=k, case=case, NAF_rep=True,
                                   family=family)

        # compute e(Q2, P2) = e(Q,P)^a (case 1) or e(Q2, P2) = e(Q,P)^ab (case 2)
        pairing_value_naf2 = ate_i(Q2, P2, P2_prec, c_vec, F, length_miller, U, W, k=k, case=case, NAF_rep=True,
                                   family=family) ** (randint_Q * randint_P)

        # verify: e(Q1, P1) = e(Q2, P2) <=> e([a]*Q,P) = e(Q,P)^a (case1)
        # verify: e(Q1, P1) = e(Q2, P2) <=> e([a]*Q,[b]*P) = e(Q,P)^ab (case2)
        print('bilinearity NAF test:',
              pairing_value_naf1 == pairing_value_naf2) if pairing_value_naf1 != 1 and pairing_value_naf2 != 1 \
            else print('review code')

    return None
