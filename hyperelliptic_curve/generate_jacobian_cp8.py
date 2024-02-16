from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ

from pairing_computation import compute_twisted_ate
from verification_operations import test_bilinearity_Twisted_Ate
from _utils import generate_curve_eq


def generate_jacobian_cp8(u, X, Y, lx, ly, l, k=8):
    """
    :param u: defines the length of the Miller loop
    :param X: used to construct the prime p as: p = X^2 + 2Y^2
    :param Y: used to construct the prime p as: p = X^2 + 2Y^2
    :param lx: used to construct the prime p as: p = X^2 + 2Y^2
    :param ly: used to construct the prime p as: p = X^2 + 2Y^2
    :param l: is used to define the twist of the Hyperelliptic curve C
    :param k: embedding degree
    :param a:
    :return:
    """
    # Hyperelliptic curve + Jacobian parameters
    r = u ** 4 + 1
    # Lifts for X and Y
    X = lx * r + X
    Y = ly * r + Y
    # Base field prime p
    p = X ** 2 + 2 * Y ** 2

    U = [u, (u // 4), lx, ly]

    # Compute the order of the Jacobian with the characteristic polynomial of Frobenius
    Zx = ZZ['t']
    (t,) = Zx._first_ngens(1)  # Define the ring of polynomials with integer coefficients
    # characteristic polynomial of Frobenius chi(t)
    xt = t ** 4 + 4 * Y * t ** 3 + 8 * Y ** 2 * t ** 2 + 4 * p * Y * t + p ** 2
    # order of the Jacobian J
    n = xt(t=1)
    # cofactor of the Jacobian J
    h = n // r
    print('Cocks-Pinch (CP) curve with embedding degree k = 8\n')
    print('---Instantiation---')
    print("u = {:#x} {} bits".format(u, u.nbits()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("h = {:#x} {} bits".format(h, h.nbits()))
    print("n = #E(Fp) = {:#x} {} bits".format(n, n.nbits()))
    print("lx = {:#x} {} bits".format(lx, lx.nbits()))

    # Construct the prime field Fp
    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']
    (x,) = Fpx._first_ngens(1)  # Fpx: ring of polynomials in x, with coefficients in Fp

    # Hyperelliptic curve C
    a = generate_curve_eq(p, n)
    C = HyperellipticCurve(x ** 5 + a * x)  # Set the equation of the hyperelliptic curve C
    # Jacobian of C
    J = C.jacobian()  # J is the Jacobian of the curve C over Fp
    F = [0, a, 0, 0, 0]

    print('\n')
    print('Hyperelleiptic curve eq. and Jacobian over the base field Fp')
    print(C)
    print(J)

    # Construct the degree 8 twist Ct of C
    # The curve Ct and its Jacobian are defined over Fp
    Ct = HyperellipticCurve(x ** 5 + a * l * x)
    Jt = Ct.jacobian()

    print('\n')
    print('Degree 8 twist of hyperelleiptic curve eq. and Jacobian over the base field Fp')
    print(Ct)
    print(Jt)

    # Construct the degree 8 extension Fp8 of the prime field Fp
    Fpw = Fp['w']
    (w,) = Fpw._first_ngens(1)  # Fpw: polynomial ring in w with coefficients in Fp
    b = 1
    while not (w ** 8 + b).is_irreducible():
        b = b + 1
    print("Fp8 = Fp[w]/(w^8+{})".format(b))  # Fp8 = Fp[w]/(w^8 + b)
    # Fp8: degree 8 extension of Fp.
    Fp8 = Fp.extension(w ** 8 + b, names=('w',))
    (w,) = Fp8._first_ngens(1)

    # Constructing the curve and Jacobian over the extension field Fp**8
    Fp8x = Fp8['x']
    (x,) = Fp8x._first_ngens(1)  # Fp8x: ring of polynomials in x, with coefficients in Fp8
    C8 = HyperellipticCurve(Fp8x([0, a, 0, 0, 0, 1]))
    J8 = C8.jacobian()

    print('\n')
    print('Hyperelleiptic curve eq. and Jacobian over the exxtension field Fp^8')
    print(C8)
    print(J8)

    # Compute the order of the Jacobian over Fp8
    Res = (t ** 8 - 1).resultant(xt)  # Compute the resultant  of the polynomials t^8 - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    # n_ = h_ * r  # Compute the order of the Jacobian J8

    curves = [C, Ct, C8]
    jacobians = [J, Jt, J8]
    fields = [Fp, Fp8]

    # Find c in Fp8 s.t. c^8 = l => c is an 8^th root of l
    g = x ** 8 - l
    roots = g.roots()
    c = roots[7]
    c = c[0]
    c_vec = []
    c_vec.append(c)
    for i in range(1, 13):
        c_vec.append(1 / (c ** i))

    print('\n')
    compute_twisted_ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, u, k=k)
    test_bilinearity_Twisted_Ate(curves, jacobians, fields, c_vec, F, U, p, r, h, h_, u, k=k)

    return None


def generate_jacobian():
    # Paremeter initialization
    # u will define the length of the Miller loop
    u = ZZ(0xffc00020fffffffc)
    assert u == 2 ** 64 - 2 ** 54 + 2 ** 37 + 2 ** 32 - 2 ** 2
    # X, Y, lx, ly are used to construct the prime p as: p = X^2 + 2Y^2
    X = ZZ(0x7fa0182f67431e596adfdc83eb3fe4757900039bffffffd8)
    Y = ZZ(0x3fc0181ce78635d69275eeda614d2265ac6f42ec69a058128c8ff985c000002d)
    lx = ZZ(0xa031)
    assert lx == 2 ** 15 + 2 ** 13 + 2 ** 5 + 2 ** 4 + 1
    ly = ZZ(1)
    # l: is used to define the twist of the Hyperelliptic curve C
    # l = 0x21272a193842552162d1c40c5258df154c31bc12353118d7e6940aa62821cbbd442344f6878da532e48272bcb2daa5a6313a0e0bbe9dabc6b450259f8f3d210aa750d39a
    l = 3
    generate_jacobian_cp8(u, X, Y, lx, ly, l)
