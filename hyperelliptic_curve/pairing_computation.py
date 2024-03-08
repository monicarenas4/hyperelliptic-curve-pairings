from write_number_operations import write_number_operations
from jacobian_operations import JC_random_element, HEC_random_point, new_coordinates, precomputation_general_div, \
    precomputation_degenerate_div
from pairing_types import twisted_ate_cp8, ate_i

file_name = 'results/number_of_operations.txt'


def compute_twisted_ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, length_miller, k: int):
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
    :param k: embedding degree
    :return:
    """
    C, Ct, Ck = curves[0], curves[1], curves[2]
    J, Jt, Jk = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fpk = fields[0], fields[1]

    # Choose random divisor P in Jacobian J(Fp)
    P = JC_random_element(C)
    P = h * P  # Force P to have order r
    # Convert P to Fan et al. coordinate system P = [U11, U10, V11, V10, Z1, Z2, z1, z2]
    # https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf
    P = new_coordinates(P)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case1: Q is a degenerate divisor: Q = [x - xQ, yQ] in Jt(Fp^s)
            Q = HEC_random_point(Ct)
            Q_prec, mult_pre, sq_pre = precomputation_degenerate_div(Q) # precomputation: values required in line evaluation
        else:
            # case2: Q is a general divisor in Mumford representation: Q = [x^2 + u21*x + u20, v21*x + v20] in Jt(Fp^s)
            Q = JC_random_element(Ct)
            Q = h_ * Q
            Q_prec, mult_pre, sq_pre = precomputation_general_div(Q) # precomputation: values required in line evaluation

        write_number_operations(file_name, embedding_degree=k, function='twisted_ate',
                                case=case, mult_pre=mult_pre, sq_pre=sq_pre, total=(mult_pre + sq_pre))

        print('----------------\nCASE: %s\n----------------' % case)
        # twisted ate pairing computation using the binary representation of length_miller
        pairing_value = twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, k=k, case=case)
        # twisted ate pairing computation using the NAF representation of length_miller
        pairing_value_naf = twisted_ate_cp8(P, Q, Q_prec, c_vec, F, length_miller, U, Fp, k=k, case=case, NAF_rep=True)
        print('pairing value = ', pairing_value)
        print('pairing value NAF = ', pairing_value_naf)

    return None


def compute_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, length_miller, k: int, family='KT16'):
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
    :param k: embedding degree
    :param family: polynomial family; used to distinguish between different constructions
    :return:
    """
    C, Ct, Ck = curves[0], curves[1], curves[2]
    J, Jt, Jk = jacobians[0], jacobians[1], jacobians[2]
    Fp, Fps, Fpk = fields[0], fields[1], fields[2]

    # Choose random divisor Q in Jacobian Jt(Fp^s)
    Q = JC_random_element(Ct)
    Q = h_ * Q  # Force Q to have order r
    # Convert Q to Fan et al. coordinate system Q = [U21, U20, V21, V20, Z1, Z2, z1, z2]
    # https://cacr.uwaterloo.ca/techreports/2008/cacr2008-03.pdf
    Q = new_coordinates(Q)

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case1: P is a degenerate divisor: P = [x - xP, yP] in J(Fp)
            P = HEC_random_point(C)
            P_prec, mult_pre, sq_pre = precomputation_degenerate_div(P) # precomputation: values required in line evaluation
        else:
            # case2: P is a general divisor in Mumford representation: P = [x^2 + u11*x + u10, v11*x + v10] in J(Fp)
            P = JC_random_element(C)
            P = h * P
            P_prec, mult_pre, sq_pre = precomputation_general_div(P) # precomputation: values required in line evaluation

        write_number_operations(file_name, embedding_degree=k, function='ate_i',
                                case=case, mult_pre=mult_pre, sq_pre=sq_pre, total=mult_pre + sq_pre)

        print('----------------\nCASE: %s\n----------------' % case)
        # ate_i pairing computation using the binary representation of length_miller
        pairing_value = ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, k=k, case=case, family=family)
        # ate_i pairing computation using the NAF representation of length_miller
        pairing_value_naf = ate_i(Q, P, P_prec, c_vec, F, length_miller, U, W, k=k,
                                  case=case, NAF_rep=True, family=family)
        print('pairing value = ', pairing_value)
        print('pairing value NAF = ', pairing_value_naf)

    return None
