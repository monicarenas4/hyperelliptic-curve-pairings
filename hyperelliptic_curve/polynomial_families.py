from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ


def polynomial_family_k16(u):
    """
    :param u: seed
    :return: r, p, X, Y
    """
    R = QQ['x']
    (x,) = R._first_ngens(1)
    rx = x ** 8 + 1
    Xx = (x ** 7 - x ** 6) / 2
    Yx = -(x ** 5 + x ** 4 + x + 1) / 4
    px = Xx ** 2 + 2 * Yx ** 2

    # Hyperelliptic curve + Jacobian parameters
    r = ZZ(rx(u) // 2)
    p = ZZ(px(u))
    X = ZZ(Xx(u))
    Y = ZZ(Yx(u))

    return r, p, X, Y


def polynomial_family_new_k16(u):
    """
    :param u: seed
    :return: r, p, X, Y
    """
    R = QQ['x']
    (x,) = R._first_ngens(1)
    rx = x ** 8 + 1
    Xx = (2 * x ** 8 + x ** 7 - x ** 6 + 2) / 2
    Yx = -(x ** 5 + x ** 4 + x + 1) / 4
    px = Xx ** 2 + 2 * Yx ** 2

    # Hyperelliptic curve + Jacobian parameters
    r = ZZ(rx(u) // 2)
    p = ZZ(px(u))
    X = ZZ(Xx(u))
    Y = ZZ(Yx(u))

    return r, p, X, Y
