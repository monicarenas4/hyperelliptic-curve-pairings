from sage.all_cmdline import *

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF

from pairings import miller_loop_tate_pairing, final_exponentiation
from checks import tate_pairing_bilinearity_check
from random import randint


def generate_curve(u: int, a: int = 0, b: int = -3):
    """
    :param u: curve seed
    :param a: coefficient a
    :param b: coefficient b
    :return: the curves (P,Q) and the integers (p,r)
    """
    R = QQ['x']
    (x,) = R._first_ngens(1)

    rx = x ** 4 - x ** 2 + 1
    tx = x + 1  # trace polynomial
    px = (x ** 6 - 2 * x ** 5 + 2 * x ** 3 + x + 1) / 3

    r = ZZ(rx(u))
    t = ZZ(tx(u))
    p = ZZ(px(u))
    n = p - t + 1  # order of the curve. #E = p - t + 1 = h*r
    h = n / r  # co-factor

    Fp = GF(p, proof=False)
    E = EllipticCurve([Fp(a), Fp(b)])
    P = E.random_element()
    P = h * P

    Fpw = Fp['w']
    (w,) = Fpw._first_ngens(1)

    # print("Fp12 = Fp[w]/(w^12 + w^6 + 2)")
    Fp12 = Fp.extension(w ** 12 + w ** 6 + 2, names=('w',))
    (w,) = Fp12._first_ngens(1)
    E12 = EllipticCurve([Fp12(a), Fp12(b)])
    n12 = E12.order()

    Q = E12.random_element()
    h12 = n12 / r ** 2
    Q = h12 * Q

    return P, Q, p, r


def tate_pairing(u):
    """
    :param u: curve seed
    :return: m, pairing value
    """
    P, Q, p, r = generate_curve(u)
    miller_function = miller_loop_tate_pairing(P, Q, r)
    pairing_value = final_exponentiation(miller_function, p, r)

    # Bilinearity check
    tate_pairing_bilinearity_check(P, Q, pairing_value, p, r)

    return miller_function, pairing_value


u = -ZZ(2 ** 63 + 2 ** 62 + 2 ** 60 + 2 ** 57 + 2 ** 48 + 2 ** 16)
tate_pairing(u)
