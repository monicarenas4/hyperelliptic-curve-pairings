from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from jacobian_operations import JC_random_element
from pairing_computation import compute_ate_i
from verification_operations import test_bilinearity_Ate_i
from _utils import w_powers_p, w_p_i, frobenius_power, generate_curve_eq
from polynomial_families import polynomial_family_k24


def generate_jacobian_k24(u, seed, k=24, family='k24'):
    """
    :param u: defines the length of the Miller loop
    :param k: embedding degree
    :param a:
    :return:
    """
    r, p, X, Y = polynomial_family_k24(u)

    print('p mod 3 = ', p % 3)

    # Compute the order of the Jacobian with the characteristic polynomial of Frobenius
    Zx = ZZ['t']
    (t,) = Zx._first_ngens(1)  # Define the ring of polynomials with integer coefficients
    # characteristic polynomial of Frobenius chi(t)
    xt = t ** 4 + 2 * (X ** 2 - 2 * Y ** 2) * t ** 2 + p ** 2
    # order of the Jacobian J
    n = xt(t=1)
    # cofactor of the Jacobian J
    h = n // r

    # Compute the order of the Jacobian over Fp8
    Res = (t ** 24 - 1).resultant(xt)  # Compute the resultant  of the polynomials t^8 - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    n_ = h_ * r  # Compute the order of the Jacobian J8

    print('---Instantiation---')
    print("u = {:#x} {} bits".format(u, u.nbits()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("h = {:#x} {} bits".format(h, h.nbits()))
    print("n = #E(Fp) = {:#x} {} bits".format(n, n.nbits()))

    # Construct the prime field Fp
    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']
    (x,) = Fpx._first_ngens(1)  # Fpx: ring of polynomials in x, with coefficients in Fp

    a = generate_curve_eq(p, n)
    print('a = ', a)
    U = [u, u - 1, 0, 0]
    F = [0, a, 0, 0, 0]

    # Hyperelliptic curve C
    C = HyperellipticCurve(x ** 5 + a * x)  # Set the equation of the hyperelliptic curve C
    # Jacobian of C
    J = C.jacobian()  # J is the Jacobian of the curve C over Fp

    Fpz = Fp['z']
    (z,) = Fpz._first_ngens(1)  # Fpw: polynomial ring in w with coefficients in Fp
    b = 1
    if seed == "u1":
        while not (z ** 3 + b + z).is_irreducible():
            b = b + 1
        print("Fp3 = Fp[z]/(z^3+ {})".format(b))  # Fp8 = Fp[w]/(w^8 + b)
        Fp3 = Fp.extension(z ** 3 + b + z, names=('z',))
    else:
        while not (z ** 3 + b).is_irreducible():
            b = b + 1
        print("Fp3 = Fp[z]/(z^3+ {})".format(b))  # Fp8 = Fp[w]/(w^8 + b)
        Fp3 = Fp.extension(z ** 3 + b, names=('z',))

    (z,) = Fp3._first_ngens(1)
    Fpw = Fp3['w']
    (w,) = Fpw._first_ngens(1)

    if seed == "u1":
        i, j = 9, 0
    else:
        i, j = 14, 0

    d = 8
    s = k // d
    print('p = ', p)
    while True:
        fX = i * z + j
        FX = w ** (k // s) - fX
        print(FX.is_irreducible())
        if FX.is_irreducible():
            l = fX
            Ct = HyperellipticCurve(x ** 5 + a * l * x)
            Jt = Ct.jacobian()

            P = JC_random_element(Ct)
            P1 = n_ * P
            P2 = h_ * P
            print(h_ * P)
            if P2[0] != 1 and P1[0] == 1:
                break
        i = i + 1
        # j = j + 1

    print('lambda = ', l)

    Fpw = Fp3['w']
    (w,) = Fpw._first_ngens(1)
    Fq8 = Fp3.extension(w ** 8 - fX, names=('w',))
    (w,) = Fq8._first_ngens(1)
    Fq8x = Fq8['x']
    (x,) = Fq8x._first_ngens(1)

    C24 = HyperellipticCurve(Fq8x([0, a, 0, 0, 0, 1]))
    J24 = C24.jacobian()

    curves = [C, Ct, C24]
    jacobians = [J, Jt, J24]
    fields = [Fp, Fp3, Fq8]

    c = Fq8.gen(0)
    # g = x ** 8 - l
    # roots = g.roots()
    # c = roots[7]
    # c = c[0]
    print('c = ', c)
    c_vec = []
    c_vec.append(c)
    for i in range(1, 13):
        c_vec.append(1 / (c ** i))

    f = Fq8.random_element()
    fp = (f ** (p ** 19))

    W = w_powers_p(w, p, k)

    fn = frobenius_power(f, k, W, 19)
    print(fp == fn)

    compute_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u, k=k, family=family)
    test_bilinearity_Ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u, k=k, family=family)

    return 0
