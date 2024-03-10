from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ

from pairing_computation import compute_ate_i
from verification_operations import test_bilinearity_ate_i
from _utils import w_powers_p
from _utils import generate_curve_eq
from polynomial_families import polynomial_family_KT16, polynomial_family_New16


def generate_jacobian_k16(u, k=16, family='KT16'):
    """
    :param u: seed to evaluate the polynomial family
    :param k: the embedding degree k = 16
    :param family: name of the family, to distinguish between KT16 and New16
    :return:
    """
    # Polynomial family for k = 16
    if family == 'KT16':
        r, p, X, Y = polynomial_family_KT16(u) # Evaluate polynomial family KT16
        U = [u, u - 1] # Values required for the final exponentation
    else:
        r, p, X, Y = polynomial_family_New16(u) # Evaluate polynomial family New16
        U = [u, u + 1] # Values required for the final exponentation

    # Compute the order of the Jacobian with the characteristic polynomial of Frobenius
    Zx = ZZ['t']
    (t,) = Zx._first_ngens(1)  # Define the ring of polynomials with integer coefficients
    # characteristic polynomial of Frobenius chi(t)
    xt = t ** 4 + 4 * Y * t ** 3 + 8 * Y ** 2 * t ** 2 + 4 * p * Y * t + p ** 2
    # order of the Jacobian J
    n = xt(t=1)
    # cofactor of the Jacobian J
    h = n // r

    print('---Instantiation---')
    print("u = {:#x} {} bits".format(u, u.nbits()))
    print("p = {:#x} {} bits".format(p, p.nbits()))
    print("r = {:#x} {} bits".format(r, r.nbits()))
    print("h = {:#x} {} bits".format(h, h.nbits()))
    print("n = #E(Fp) = {:#x} {} bits".format(n, n.nbits()))

    # Determine the curve coefficient a
    if family == 'KT16':
        a = generate_curve_eq(p, n)
    else:
        a = generate_curve_eq(p, n)

    # Construct the prime field Fp
    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']
    (x,) = Fpx._first_ngens(1)  # Fpx: ring of polynomials in x, with coefficients in Fp

    # Hyperelliptic curve C
    C = HyperellipticCurve(x ** 5 + a * x)  # Set the equation of the hyperelliptic curve C
    # Jacobian of C
    J = C.jacobian()  # J is the Jacobian of the curve C over Fp
    # F vector contains coefficients of the polynomial f(x) s.t. y^2 = f(x)
    # this is needed in the divisor doubling and addition
    # in our case, this is not used in doubling and addition because most coefficients are 0
    F = [0, a, 0, 0, 0]

    print('\n')
    print('Hyperelliptic curve eq. and Jacobian over the base field Fp')
    print(C)
    print(J)

    # Construct the degree 2 extension Fq = Fp2 of the prime field Fp
    Fpz = Fp['z']
    (z,) = Fpz._first_ngens(1)  # Fpz: polynomial ring in z with coefficients in Fp
    b = 1
    while not (z ** 2 + b).is_irreducible():
        b = b + 1
    print("Fp2 = Fp[z]/(z^2+{})".format(b))  # Fq = Fp[z]/(z^2 + b)

    # Fq = Fp2: degree 2 extension of Fp.
    Fp2 = Fp.extension(z ** 2 + b, names=('z',))
    (z,) = Fp2._first_ngens(1)

    # Define the twist parameter l to define the degree 8 twist Ct: y^2 = x^5 + a*l*x of C
    # TODO: modify this to work for non-fixed values of l
    if family == 'KT16':
        l = z
    else:
        l = 5*z

    # Construct the degree 8 twist Ct of C. Ct is defined over Fp
    # Construct the Jacobian Jt of the twist Ct. Jt is defined over Fp
    Ct = HyperellipticCurve(x ** 5 + a * l * x)
    Jt = Ct.jacobian()

    print('\n')
    print('Degree 8 twist of hyperelliptic curve eq. and Jacobian over the extension field Fq = Fp^2')
    print(Ct)
    print(Jt)

    # Construct the degree 8 extension Fq8 of the field Fq = Fp^2
    Fpw = Fp2['w']
    (w,) = Fpw._first_ngens(1)
    Fq8 = Fp2.extension(w ** 8 - l, names=('w',)) # Fq^8 = Fp^2[w]/(w^8 - l)
    (w,) = Fq8._first_ngens(1)
    Fq8x = Fq8['x']
    (x,) = Fq8x._first_ngens(1)

    # Construct the curve C16: y^2 = x^5 + a*x and Jacobian J16 over the extension field Fq^8
    C16 = HyperellipticCurve(Fq8x([0, a, 0, 0, 0, 1]))
    J16 = C16.jacobian()

    print('\n')
    print('Hyperelliptic curve eq. and Jacobian over the extension field Fq^8')
    print(C16)
    print(J16)

    # Compute the order of the Jacobian over Fq^8
    Res = (t ** k - 1).resultant(xt)  # Compute the resultant  of the polynomials t^k - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    n_ = h_ * r  # Compute the order of the Jacobian J8

    curves = [C, Ct, C16]
    jacobians = [J, Jt, J16]
    fields = [Fp, Fp2, Fq8]

    # The values in c_vec are needed to move from Fq to Fq8 in line evaluation
    c = Fq8.gen(0)
    c_vec = []
    c_vec.append(c)
    for i in range(1, 13):
        c_vec.append(1 / (c ** i))

    # powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    W = w_powers_p(c, p, k)

    print('\n')
    compute_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u, k=k, family=family)
    test_bilinearity_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u, k=k, family=family)

    return None
