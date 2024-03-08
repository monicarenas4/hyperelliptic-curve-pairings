from sage.all import Integer
from math import log2, floor

from jacobian_operations import ADD, DBL
from _utils import NAF, hamming_weight, NAf_hamming_weight, field_conversion
from write_number_operations import write_number_operations

file_name = 'results/number_of_operations.txt'


def miller_function(P, Q, Q_prec, c_vec, F, length_miller, case: str, k: int, twist: bool = False, NAF_rep=False):
    """
    :param P: divisor over Fp (for twisted ate pairing) or over Fp^s (for ate_i) in Fan et al. coordinate system
    :param Q: divisor over Fp^s (for twisted ate pairing) or over Fp (for ate_i) represented as degenerate/general divisor
    :param Q_prec: precomputation values required for line evaluation in Miller loop
    :param c_vec: powers of c, the generator of the extension field Fq^d, used in line evaluation
    :param F: coefficients of the polynomial f(x), s.t. C/Fp: y^2 = f(x)
    :param length_miller: integer that defined the length of the Miller loop. In all cases this is the seed u
    :param case: case1 => degenerate divisor or case2 => general divisor
    :param twist: distinguish between twisted ate and ate_i pairings
    :param NAF_rep: determine whether NAF or binary representation is used
    :return: fc, the Miller function
    """
    global P_neg, mult_DBL, sq_DBL, mult_ADD, sq_ADD, const_mult_line_DBL, mult_line_DBL, sq_line_DBL
    global const_mult_line_ADD, mult_line_ADD, sq_line_ADD, sk, mk, s, m, cm, mk_DBL, mult_miller_DBL

    if not NAF_rep:
        # binary representation of length_miller: vector_miller = (b_n, ..., b_0) with b_i = {0,1}
        vector_miller = Integer(length_miller).digits(2)
    else:
        # NAF representation of length_miller: vector_miller = (b_n, ..., b_0) with b_i = {-1,0,1}
        vector_miller = NAF(length_miller)
        # divisor -P in Fan et al. coordinate system
        P_neg = [P[0], P[1], -P[2], -P[3], P[4], P[5], P[6], P[7]]

    T, fc = P, 1
    number_of_DBL = length_miller.nbits() - 1 # number_of_DBL = log(length_miller) - 1 = n - 1: number of doubling steps

    for i in range(len(vector_miller) - 2, -1, -1):
        # Doubling step. T = [2]T and lc = line computation
        T, lc, const_mult_line_DBL, mult_line_DBL, sq_line_DBL, mult_DBL, sq_DBL = DBL(T, Q_prec, Q, F, c_vec,
                                                                                       case=case, twist=twist)
        # DBL update for Miller function
        fc = lc * fc ** 2  # 1_m8, 1_s8 (sparse)
        if vector_miller[i] == 1:
            # Addition step. T = T + P and lc = line computation
            T, lc, const_mult_line_ADD, mult_line_ADD, sq_line_ADD, mult_ADD, sq_ADD = ADD(P, T, Q_prec, Q, F, c_vec,
                                                                                           case=case, twist=twist)
            # ADD update for Miller function
            fc = lc * fc  # 1_m8 (sparse)
        elif vector_miller[i] == -1:
            # Addition step. T = T - P and lc = line computation (for bits b_i = -1)
            T, lc, const_mult_line_ADD, mult_line_ADD, sq_line_ADD, mult_ADD, sq_ADD = ADD(P_neg, T, Q_prec, Q, F,
                                                                                           c_vec, case=case,
                                                                                           twist=twist)
            # ADD update for Miller function
            fc = lc * fc  # 1_m8

    if not NAF_rep:
        # number of additions = hw(length_miller) - 1
        number_of_ADD = hamming_weight(vector_miller) - 1
    else:
        # number of additions = hw_NAF(length_miller) - 1
        number_of_ADD = NAf_hamming_weight(vector_miller) - 1

    m, s, cm, mk, mk_DBL, sk = field_conversion(k)

    mult_DBL = (((number_of_DBL - 1) * (mult_DBL * m + const_mult_line_DBL * cm + mult_line_DBL * m))
                + (1 * 25 * m + const_mult_line_DBL * cm + mult_line_DBL * m))
    sq_DBL = ((number_of_DBL - 1) * (sq_DBL * s + sq_line_DBL * s)) + (1 * 5 * s + sq_line_DBL * s)
    mult_ADD = number_of_ADD * (mult_ADD * m + const_mult_line_ADD * cm + mult_line_ADD * m)
    sq_ADD = number_of_ADD * (sq_ADD * s + sq_line_ADD * s)
    mult_miller_DBL = (mk_DBL * (number_of_DBL - 1))
    sq_miller_DBL = (sk * (number_of_DBL - 1))
    mult_miller_ADD = (mk_DBL * number_of_ADD)
    sq_miller_ADD = 0 * number_of_ADD

    total = mult_DBL + sq_DBL + mult_ADD + sq_ADD + mult_miller_DBL + sq_miller_DBL + mult_miller_ADD + sq_miller_ADD

    write_number_operations(file_name,
                            embedding_degree=k,
                            function='miller', case=case, twist=twist, NAF_rep=NAF_rep,
                            number_of_DBL=number_of_DBL,
                            number_of_ADD=number_of_ADD,
                            mult_DBL=mult_DBL, sq_DBL=sq_DBL,
                            mult_ADD=mult_ADD, sq_ADD=sq_ADD,
                            mult_miller_DBL=mult_miller_DBL,
                            sq_miller_DBL=sq_miller_DBL,
                            mult_miller_ADD=mult_miller_ADD,
                            sq_miller_ADD=sq_miller_ADD, total=total)

    return fc
