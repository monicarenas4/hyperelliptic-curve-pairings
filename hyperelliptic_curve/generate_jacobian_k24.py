from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

from jacobian_operations import JC_random_element
from pairing_computation import compute_ate_i
from verification_operations import test_bilinearity_ate_i
from _utils import w_powers_p, w_p_i, frobenius_power, generate_curve_eq
from polynomial_families import polynomial_family_New24


def generate_jacobian_k24(u, seed, k=24, family='New24'):
    """
    :param u: seed to evaluate the polynomial family
    :param seed: name of the seed to distinguish between different instantiations
    :param k: the embedding degree k = 24
    :param family: name of the family
    :return:
    """

    # Evaluate polynomial family New24
    r, p, X, Y = polynomial_family_New24(u)

    # Compute the order of the Jacobian with the characteristic polynomial of Frobenius
    Zx = ZZ['t']
    (t,) = Zx._first_ngens(1)  # Define the ring of polynomials with integer coefficients
    # characteristic polynomial of Frobenius chi(t)
    xt = t ** 4 + 2 * (X ** 2 - 2 * Y ** 2) * t ** 2 + p ** 2
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

    # Construct the prime field Fp
    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']
    (x,) = Fpx._first_ngens(1)  # Fpx: ring of polynomials in x, with coefficients in Fp

    # Determine the curve coefficient a
    a = generate_curve_eq(p, n)
    # Values required for the final exponentation
    U = [u, u - 1]

    # F vector contains coefficients of the polynomial f(x) s.t. y^2 = f(x)
    # this is needed in the divisor doubling and addition
    # in our case, this is not used in doubling and addition because most coefficients are 0
    F = [0, a, 0, 0, 0]

    # Hyperelliptic curve C
    C = HyperellipticCurve(x ** 5 + a * x)  # Set the equation of the hyperelliptic curve C
    # Jacobian of C
    J = C.jacobian()  # J is the Jacobian of the curve C over Fp

    print('\n')
    print('Hyperelliptic curve eq. and Jacobian over the base field Fp')
    print(C)
    print(J)

    # Construct the degree 3 extension Fq = Fp^3 of the prime field Fp
    Fpz = Fp['z']
    (z,) = Fpz._first_ngens(1)  # Fpz: polynomial ring in z with coefficients in Fp
    b = 1
    if seed == "u1":
        while not (z ** 3 + z + b).is_irreducible():
            b = b + 1
        # in this case, p != 1 mod 3 so the degree 3 ext. field is: # Fp^3 = Fp[z]/(z^3 + z + b)
        print("Fp3 = Fp[z]/(z^3 + z + {})".format(b))
        Fp3 = Fp.extension(z ** 3 + b + z, names=('z',))
    else:
        while not (z ** 3 + b).is_irreducible():
            b = b + 1
        # in this case, p = 1 mod 3 so the degree 3 ext. field is: # Fp^3 = Fp[z]/(z^3 + b)
        print("Fp3 = Fp[z]/(z^3 + {})".format(b))
        Fp3 = Fp.extension(z ** 3 + b, names=('z',))
    (z,) = Fp3._first_ngens(1)

    # Define the twist parameter l to define the degree 8 twist Ct: y^2 = x^5 + a*l*x of C
    # TODO: modify this to work for non-fixed values of l
    if seed == "u1":
        l = 9*z
    else:
        l = 14*z

    # Construct the degree 8 twist Ct of C. Ct is defined over Fp
    # Construct the Jacobian Jt of the twist Ct. Jt is defined over Fp
    Ct = HyperellipticCurve(x ** 5 + a * l * x)
    Jt = Ct.jacobian()

    print('\n')
    print('Degree 8 twist of hyperelliptic curve eq. and Jacobian over the extension field Fq = Fp^3')
    print(Ct)
    print(Jt)

    # Construct the degree 8 extension Fq8 of the field Fq = Fp^3
    Fpw = Fp3['w']
    (w,) = Fpw._first_ngens(1)
    Fq8 = Fp3.extension(w ** 8 - l, names=('w',)) # Fq^8 = Fp^3[w]/(w^8 - l)
    (w,) = Fq8._first_ngens(1)
    Fq8x = Fq8['x']
    (x,) = Fq8x._first_ngens(1)

    # Construct the curve C24: y^2 = x^5 + a*x and Jacobian J24 over the extension field Fq^8
    C24 = HyperellipticCurve(Fq8x([0, a, 0, 0, 0, 1]))
    J24 = C24.jacobian()

    print('\n')
    print('Hyperelliptic curve eq. and Jacobian over the extension field Fq^8')
    print(C24)
    print(J24)

    # Compute the order of the Jacobian over Fq^8
    Res = (t ** k - 1).resultant(xt)  # Compute the resultant  of the polynomials t^24 - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    n_ = h_ * r  # Compute the order of the Jacobian J8

    curves = [C, Ct, C24]
    jacobians = [J, Jt, J24]
    fields = [Fp, Fp3, Fq8]

    # The values in c_vec are needed to move from Fq to Fq8 in line evaluation
    c = Fq8.gen(0)
    c_vec = []
    c_vec.append(c)
    for i in range(1, 13):
        c_vec.append(1 / (c ** i))

    # powers of p^j of the generator w of Fq8 - needed for computing Frobenius powers for elements in Fq^8
    W = w_powers_p(w, p, k)

    compute_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u, k=k, family=family)
    test_bilinearity_ate_i(curves, jacobians, fields, c_vec, F, U, W, p, r, h, h_, u, k=k, family=family)

    return 0
