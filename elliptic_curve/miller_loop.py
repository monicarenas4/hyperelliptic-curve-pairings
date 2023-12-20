from sage.all_cmdline import *
from elliptic_curve.ecc_operations import ADD_step, DBL_step, ADD_miller_point, DBL_miller_point


def miller_loop_tate_pairing(P, Q, r: int):
    """
    :param P: point on the curve E(Fp)
    :param Q: point on the extension filed E(Fp12)
    :param r: prime number that divides the oder of the curve
    :return: miller_function
    """
    bin_r = [int(x) for x in "{0:0b}".format(r)][1:]
    R = P
    miller_function = 1

    for bit in bin_r:
        R, l, v = DBL_step(Q, R)
        miller_function = miller_function ** 2 * (l / v)
        if bit == 1:
            R, l, v = ADD_step(P, Q, R)
            miller_function = miller_function * (l / v)

    return miller_function


def miller_loop_opt_tate_pairing(P, Q, r: int):
    """
    :param P: point on the curve E(Fp)
    :param Q: point on the extension filed E(Fp12)
    :param r: prime number that divides the oder of the curve
    :return: miller_function
    """
    bin_r = [int(x) for x in "{0:0b}".format(r)][1:]
    R = P
    miller_function_num = 1
    miller_function_den = 1

    for bit in bin_r:
        R, l, v = DBL_miller_point(Q, R)
        miller_function_num = miller_function_num ** 2 * l
        miller_function_den = miller_function_den ** 2 * v
        if bit == 1:
            R, l, v = ADD_miller_point(P, Q, R)
            miller_function_num = miller_function_num * l
            miller_function_den = miller_function_den * v

    miller_function = miller_function_num / miller_function_den

    return miller_function


def miller_loop_opt_ate_pairing(Q, P, u: int):
    """
    :param P: point on the curve E(Fp)
    :param Q: point on the extension filed E(Fp12)
    :param r: prime number that divides the oder of the curve
    :return: miller_function
    """
    bin_r = [int(x) for x in "{0:0b}".format(abs(u))][1:]
    R = Q
    miller_function_num = 1
    miller_function_den = 1

    for bit in bin_r:
        R, l, v = DBL_miller_point(P, R)
        miller_function_num = miller_function_num ** 2 * l
        miller_function_den = miller_function_den ** 2 * v
        if bit == 1:
            R, l, v = ADD_miller_point(Q, P, R)
            miller_function_num = miller_function_num * l
            miller_function_den = miller_function_den * v

    if u < 0:
        miller_function = miller_function_den / miller_function_num
    else:
        miller_function = miller_function_num / miller_function_den

    return miller_function
