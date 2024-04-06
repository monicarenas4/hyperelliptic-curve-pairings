from sage.all_cmdline import *

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF


def BLS12_curve(u: int, a: int = 0, b: int = -3, type: str = None):
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
    if type == 'ate':
        xQ_prime, yQ_prime = Q[0].frobenius(), Q[1].frobenius()
        Q_prime = E12(xQ_prime, yQ_prime)
        Q = Q_prime - Q

    return P, Q, p, r

def FK12_curve(u: int):
    R = QQ['x']
    (x,) = R._first_ngens(1)

    rx = 36 * x ** 4 + 36 * x ** 3 + 18 * x ** 2 + 6 * x + 1
    tx = -6 * x ** 2 + 1  # trace polynomial
    px = 1728 * x ** 6 + 2160 * x ** 5 + 1548 * x ** 4 + 756 * x ** 3 + 240 * x ** 2 + 54 * x + 7

    r = ZZ(rx(u))
    t = ZZ(tx(u))
    p = ZZ(px(u))
    n = p - t + 1  # order of the curve. #E = p - t + 1 = h*r
    h = n / r  # co-factor

    print("u = {:#x} {} bits".format(u, u.nbits()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))

    print(is_prime(r))
    print(is_prime(p))