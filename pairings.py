from sage.all_cmdline import *
from sage.rings.integer import Integer
from ecc_operations import ADD_step, DBL_step


def tate_pairing_miller_loop(P, Q, r: int):
    """
    :param P: point on the curve E(Fp)
    :param Q: point on the extension filed E(Fp12)
    :param r: prime number that divides the oder of the curve
    :return: None
    """
    bin_r = Integer(r).digits(2)
    R = P
    miller_function = 1

    for bit in bin_r[1:]:
        R, l, v = DBL_step(Q, R)
        miller_function = miller_function ** 2 * (l / v)
        if bit == 1:
            l, v = ADD_step(P, Q, R)
            miller_function = miller_function * (l / v)

    return miller_function
