from operations import JC_random_element, HEC_random_point, precomputation_general_div
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ
from _utils import NAF
from operations import new_coordinates
from twisted_operations import Twisted_Ate_k8
from random import randint


def generate_jacobian(k: int = 8, a: int = 3):
    """
    :param k: embedding degree
    :param a:
    :return:
    """
    p = 0x63757e4ba8c6ff6428754c24f70a9d3fe49534369f934a87b3bc1ff278656337cc69cee396a8ef98ad875836188ff0f293ae2a233bd903541cf070deadb7631ff5f27ad9
    r = 0xff0060739e18d7594a978b0ab6ae4ce3dbfd52a9d00197603fffdf0000000101
    h = 0x26cad1cf6fb1762e04b1549002acb3556aa8178f23bd901d1d01f940fb055fb7ca43e8b854a30786a394a65690a583fbb88c4c850a7fcf78daf75074603484a1c06a742ea4a9d002bf9b63808aeee5759acee12b509649987d7270d3c561273221ebfbba91d5a0c2
    n = 0x26a4159b3400bfed201dba82df5c3f55f75b70984638ccda45d4079e5d3b97c8b78d2bb8fbff84726afe91c6c4112fb96ca0a1716c12a0eae299b835cd4c05623913386752579775193e447b5ebf1b530b78dc7b5bcfedfb337885eae68ea3a4b994ee7ea2a443d8c1daf95bc29d0b37b8037ae7968df83ff7c1a7a9523b6b78042eb44c677662c2

    # X = 0x9f910b5ab8f0c95906b67d8867694eae3f74da06d105ba2a05e1a459b29c00a0d109
    Y = 0x13ec07890859f0d2fdd0d79e517fb6f49886c959639a1ef72cc8fd885c000012e
    T = 0xffc00020fffffffc

    l = 7457621329808235781107228406086954129386944851397679741565717026992476628562587060421593049761741274109849681179124230276920918105899662131351644732718514068050842

    NAF_T = NAF(T)

    Zx = ZZ['t']  # Define the ring of polynomials with integer coefficients
    (t,) = Zx._first_ngens(1)
    xt = t ** 4 + 4 * Y * t ** 3 + 8 * Y ** 2 * t ** 2 + 4 * p * Y * t + p ** 2
    assert n == xt(t=1)  # Check is the order is correct

    lx, ly = 41009, 1
    U = [T, (T // 4), lx, ly]
    F = [0, 3, 0, 0, 0]

    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']  # Fix Fpx as the ring of polynomials in x, with coefficients in Fp
    (x,) = Fpx._first_ngens(1)

    # Hyperelliptic curve
    C = HyperellipticCurve(x ** 5 + 3 * x)  # Set the equation of the hyperelliptic curve C
    # Jac = C.jacobian()  # J is the Jacobian of the curve C

    P = JC_random_element(C)
    P = h * P  # Force P to have order r

    randint_a = randint(0, r - 1)
    M = randint_a * P

    # R1, R2 = P + M, 2 * P

    D1 = new_coordinates(P)
    D1_ = -P
    D1_ = new_coordinates(D1_)

    D2 = new_coordinates(M)
    D2_ = -M
    D2_ = new_coordinates(D2_)

    Ct = HyperellipticCurve(x ** 5 + 3 * l * x)
    # Jt = Ct.jacobian()

    # Constructing the curve and Jacobian over the extension field Fp**8
    # print("Fp12 = Fp[w]/(w^12 + w^6 + 2)")

    # Fpw: polynomial ring in w with coefficients in Fp
    Fpw = Fp['w']
    (w,) = Fpw._first_ngens(1)
    Fp8 = Fp.extension(w ** k + a, names=('w',))  # Fp8: degree 8 extension of Fp.
    (w,) = Fp8._first_ngens(1)
    Fp8x = Fp8['x']
    (x,) = Fp8x._first_ngens(1)
    C8 = HyperellipticCurve(Fp8x([0, 3, 0, 0, 0, 1]))
    Jac8 = C8.jacobian()
    Res = (t ** 8 - 1).resultant(xt)  # Compute the resultant  of the polynomials t^8 - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    n_ = h_ * r  # Compute the order of the Jacobian J8

    g = x ** 8 - l
    roots = g.roots()
    c = roots[7]
    c = c[0]
    L = [1, 1, 1, 1]

    cases = ['case1', 'case2']

    for case in cases:
        if case == 'case1':
            # case 1 -> Degenerate Divisor
            randint_b = 1
            Q = HEC_random_point(Ct)
            xQ, yQ = Q[0], Q[1]

            Q = Ct([xQ, yQ])

            xQ, yQ = xQ / (c ** 2), yQ / (c ** 5)

            Q = C8([xQ, yQ])

            # Q = [-xQ, yQ, xQ**2, -xQ**3]
            Q1, Q2 = [-xQ, yQ], [-xQ, yQ]
            Q1_vec, Q2_vec = [xQ ** 2, -xQ ** 3], [xQ ** 2, -xQ ** 3]
        else:
            # case 2 -> General divisor
            Q1 = HEC_random_point(Ct)
            Q2 = HEC_random_point(Ct)
            xQ1, yQ1 = Q1[0], Q1[1]
            xQ2, yQ2 = Q2[0], Q2[1]

            xQ1 = xQ1 / (c ** 2)
            yQ1 = yQ1 / (c ** 5)
            xQ2 = xQ2 / (c ** 2)
            yQ2 = yQ2 / (c ** 5)

            Q1 = C8([xQ1, yQ1])
            Q2 = C8([xQ2, yQ2])

            u21, u20 = -(xQ1 + xQ2), (xQ1 * xQ2)

            v21 = (yQ2 - yQ1) / (xQ2 - xQ1)
            v20 = ((xQ2 - xQ1) * yQ1 - (yQ2 - yQ1) * xQ1) / (xQ2 - xQ1)

            ux, vx = (x ** 2 + u21 * x + u20), (v21 * x + v20)

            Q1 = Jac8([ux, vx])
            Q1 = h_ * Q1

            Q1_vec = precomputation_general_div(Q1)

            randint_b = randint(0, r - 1)
            Q2 = randint_b * Q1
            Q2_vec = precomputation_general_div(Q2)

        t1 = Twisted_Ate_k8(D=D2, Dn=D2_, Q_vec=Q2_vec, Q=Q2, F=F, s=NAF_T, U=U, K=Fp, L=L, case=case)
        t2 = Twisted_Ate_k8(D=D1, Dn=D1_, Q_vec=Q1_vec, Q=Q1, F=F, s=NAF_T, U=U, K=Fp, L=L, case=case) ** (
                    randint_a * randint_b)
        print(case)
        print(t1)
        print(t2)
        print(t1 == t2)

    return U


generate_jacobian()
