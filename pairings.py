from sage.all_cmdline import *
from ecc_operations import ADD_step, DBL_step


def tate_pairing_miller_loop(P, Q, r: int):
    """
    P: point on the curve E(Fp)
    Q: point on the extension filed E(Fp12)
    r: prime that divides the oder of the curve
    """
    r = int('{0:b}'.format(r))
    R = P
    miller_function = 1

    for bin_i in r[1:]:
        R, l, v = DBL_step(P, R, Q)
        miller_function = miller_function ** 2 * (l / v)
        if bin_i == 1:
            l, v = ADD_step(P, R, Q)
            miller_function = miller_function * (l / v)

    return miller_function
