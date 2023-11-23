from sage.all_cmdline import *
from sage.rings.integer import Integer
from ecc_operations import ADD_step, DBL_step


def tate_pairing_miller_loop(P, Q, r: int):
    """
    :param P: point on the curve E(Fp)
    :param Q: point on the extension filed E(Fp12)
    :param r: prime number that divides the oder of the curve
    :return: miller_function
    """
    bin_r = Integer(r).digits(2)
    bin_r.reverse()
    print(len(bin_r))
    R = P
    miller_function = 1
    count = 0

    for bit in bin_r[1:]:
        count += 1
        print(count)
        R, l, v = DBL_step(Q, R)
        print("R\t", R) if count == 254 else None
        miller_function = miller_function ** 2 * (l / v)
        if bit == 1:
            R, l, v = ADD_step(P, Q, R)
            miller_function = miller_function * (l / v)
    return miller_function


def final_exponentiation(miller_fn, p, r, k: int = 12):
    """
    :param miller_fn: miller function
    :param p: prime number
    :param r: prime number
    :param k: embedding factor
    :return: exponentiation
    """
    exponentiation = miller_fn ** ((p ** k - 1) / r)

    return exponentiation
