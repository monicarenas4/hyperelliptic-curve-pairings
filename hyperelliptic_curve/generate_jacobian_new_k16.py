from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from pairing_computation import test_ate_i
from verification_operations import test_bilinearity_Ate_i
from _utils import w_powers_p, w_p_i, frobenius_power
from _utils import generate_curve_eq

def generate_jacobian_new_k16(u, k=16, a=7):
    """
    :param u: defines the length of the Miller loop
    :param k: embedding degree
    :param a:
    :return:
    """
    # Polynomial family for k = 16
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
    
    print("p = {:#x} {}-bits".format(p, p.nbits()))
    print("r = {:#x} {}-bits".format(r, r.nbits()))
    
    U = [u, u + 1, 0, 0]

    # Compute the order of the Jacobian with the characteristic polynomial of Frobenius
    Zx = ZZ['t']
    (t,) = Zx._first_ngens(1)  # Define the ring of polynomials with integer coefficients
    # characteristic polynomial of Frobenius chi(t)
    xt = t ** 4 + 4 * Y * t ** 3 + 8 * Y ** 2 * t ** 2 + 4 * p * Y * t + p ** 2
    # order of the Jacobian J
    n = xt(t=1)
    # cofactor of the Jacobian J
    h = n // r

    # Construct the prime field Fp
    a = generate_curve_eq(p, n)
    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']
    (x,) = Fpx._first_ngens(1)  # Fpx: ring of polynomials in x, with coefficients in Fp

    a = 29
    # Hyperelliptic curve C
    C = HyperellipticCurve(x ** 5 + a * x)  # Set the equation of the hyperelliptic curve C
    # Jacobian of C
    J = C.jacobian()  # J is the Jacobian of the curve C over Fp
        
    F = [0, a, 0, 0, 0]

    Fpz = Fp['z']
    (z,) = Fpz._first_ngens(1)  # Fpw: polynomial ring in w with coefficients in Fp
    b = 1
    while not (z ** 2 + b).is_irreducible():
        b = b + 1
    print("Fp2 = Fp[z]/(z^2+{})".format(b))  # Fp8 = Fp[w]/(w^8 + b)
    # Fp8: degree 8 extension of Fp.
    Fp2 = Fp.extension(z ** 2 + b, names=('z',))
    (z,) = Fp2._first_ngens(1)

    l = z
    Ct = HyperellipticCurve(x ** 5 + a * l * x)
    Jt = Ct.jacobian()

    Fpw = Fp2['w']
    (w,) = Fpw._first_ngens(1)
    Fq8 = Fp2.extension(w ** 8 - z, names=('w',))
    (w,) = Fq8._first_ngens(1)
    Fq8x = Fq8['x']
    (x,) = Fq8x._first_ngens(1)

    C16 = HyperellipticCurve(Fq8x([0, a, 0, 0, 0, 1]))
    J16 = C16.jacobian()

    # Compute the order of the Jacobian over Fp8
    Res = (t ** 16 - 1).resultant(xt)  # Compute the resultant  of the polynomials t^8 - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    n_ = h_ * r  # Compute the order of the Jacobian J8

    curves = [C, Ct, C16]
    jacobians = [J, Jt, J16]
    fields = [Fp, Fp2, Fq8]

    c = Fq8.gen(0)
    c_vec = []
    for i in range(1, k):
        c_vec.append(c ** i)

    W = w_powers_p(w, p, k)
    
#    test_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u)
    test_bilinearity_Ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u)

    return None
