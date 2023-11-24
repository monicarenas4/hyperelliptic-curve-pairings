from sage.all_cmdline import *

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF

from pairings import miller_loop_tate_pairing, final_exponentiation
from random import randint


def tate_pairing_bilinearity_check(P, Q, pairing_value, p: int, r: int):
    """
    :param P:
    :param Q:
    :param p:
    :param r:
    :return:
    """
    a, b = randint(0, r - 1), randint(0, r - 1)
    pairing_value = pairing_value ** (a * b)

    # pairing of P' and Q'
    P_prime, Q_prime = a * P, b * Q
    miller_function_prime = miller_loop_tate_pairing(P_prime, Q_prime, r)
    pairing_value_prime = final_exponentiation(miller_function_prime, p, r)

    print("Bilinearity Check:", (pairing_value == pairing_value_prime))

    return None
