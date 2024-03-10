import os
import time
from sage.all import Integer
from jacobian_operations import JC_random_element
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve


def make_folder(folder_name: str):
    """
    :param folder_name: folder name
    :return:
    """
    os.makedirs(folder_name) if not os.path.exists(folder_name) else None
    return


def NAF(x: int):
    """
    :param x: integer x
    :return: list corresponding to the NAF represnetation of x
    """
    t0 = time.time()

    naf_x = []
    xx = Integer(x)
    assert x >= 0
    while xx > 0:
        rr = xx % 4
        if rr == 3:
            rr = -1
        else:
            rr = rr % 2
        naf_x.append(rr)
        xx -= rr
        xx, rr = xx.quo_rem(2)
        assert rr == 0
    assert x == sum([r * 2 ** i for i, r in enumerate(naf_x)])

    return naf_x


def hamming_weight(bit_x: list) -> int:
    """
    :param bit_x: binary representartion of a positive integer x
    :return: Hamming weight of x
    """
    count = 0
    for i in range(len(bit_x)):
        if bit_x[i] == 1:
            count = count + 1

    return count


def NAf_hamming_weight(naf_x: list) -> int:
    """
    :param naf_x: NAF representation of a positive integer x
    :return: NAF Hamming weight of x
    """
    count = 0
    for i in range(0, len(naf_x)):
        if (naf_x[i] == 1) or (naf_x[i] == -1):
            count = count + 1

    return count


def w_p_i(w, p, k, j):
    """
    :param w: generator of the extension field Fq^d, where Fq = Fp^s
    :param p: base field prime
    :param k: embedding degree
    :param j: a positive integer j < k
    :return: w_power a list of powers of w raised tp p^j

    the function computes the list: w_power = [w^(p^j), w^(2p^j), ..., w^((k - 1)p^j)]
    """
    w_power = []
    for i in range(0, k):
        w_ = (w ** i) ** (p ** j)
        w_power.append(w_)
    return w_power


def w_powers_p(w, p, k):
    """
    :param w: generator of the extension field Fq^d, where Fq = Fp^s
    :param p: base field prime
    :param k: embedding degree
    :return: W a list of powers of w raised tp p^j, for every j = 0, ..., k - 1

    the function computes the list: W = [w_power_1, w_power_2, ..., w_power_(k-1)], where:

    w_power_1 = [w^p, w^2p, ..., w^(k - 1)p]
    w_power_2 = [w^(p^2), w^(2p^2), ..., w^((k - 1)p^2)]
    ...
    w_power_2 = [w^(p^(k - 1)), w^(2p^(k - 1)), ..., w^((k - 1)p^(k - 1))]
    """
    W = []
    for j in range(0, k):
        w_pow = w_p_i(w, p, k, j)
        W.append(w_pow)
    return W


def frobenius_power(f, k, W, d, j):
    """
    :param f: element f in Fq^d = Fp^k
    :param k: embedding degree
    :param W: a (precomputed) list of powers of w raised tp p^j, for every j = 0, ..., k - 1
    :param d: the degree of the twist
    :param j: a positive integer j < k indicating the power of p that will be computed with Frobenius
    :return: the element f^(p^j)
    """
    s = k // d
    k_ = k // s
    f_coeff = []
    for i in range(0, k_):
        f_frob = (f[i]).frobenius(j)
        f_coeff.append(f_frob)
    f_power = f_coeff[0]
    for i in range(1, k_):
        f_power = f_power + (f_coeff[i] * W[j][i])
    return f_power


def generate_curve_eq(p, n):
    """
    :param p: base field prime
    :param n: the order of Jacobian J(Fp)
    :return: coefficient a of the curve C/Fp

    The function works only for curves of the form: C/Fp: y^2 = x^5 + a*x
    """
    Fp = GF(p, proof=False)     # Fix the prime field Fp
    Fpx = Fp['x']
    (x,) = Fpx._first_ngens(1)  # Fpx: ring of polynomials in x, with coefficients in Fp
    a = 0
    while True:
        a = a + 1
        C = HyperellipticCurve(x ** 5 + a * x)
        J = C.jacobian()
        P = JC_random_element(C)
        P = n * P
        if P[0] == 1:
            return a

def field_conversion(k: int):
    """
    :param k: embedding degree
    :return: m, s, cm, mk, mk_DBL, sk
    """

    m, s, cm, mk, mk_DBL, sk = 1, 1, 1, 1, 1, 1
    if k == 8:
        m1, s1 = 1, 1
        m, s, cm, mk, mk_DBL, sk = (m * m1), (s * s1), 1, 27, 18, 18
    elif k == 16:
        m2, s2 = 3, 2
        m, s, cm, mk, mk_DBL, sk = (m * m2), (s * s2), 2, 81, 54, 54
    elif k == 24:
        m3, s3 = 6, 5
        m, s, cm, mk, mk_DBL, sk = (m * m3), (s * s3), 3, 162, 108, 108

    return m, s, cm, mk, mk_DBL, sk