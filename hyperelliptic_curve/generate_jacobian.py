from operations import JC_random_element
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.rings.integer_ring import ZZ
from _utils import NAF


def generate_jacobian(k=8, a=3):
    p = 0x63757e4ba8c6ff6428754c24f70a9d3fe49534369f934a87b3bc1ff278656337cc69cee396a8ef98ad875836188ff0f293ae2a233bd903541cf070deadb7631ff5f27ad9
    r = 0xff0060739e18d7594a978b0ab6ae4ce3dbfd52a9d00197603fffdf0000000101
    h = 0x26cad1cf6fb1762e04b1549002acb3556aa8178f23bd901d1d01f940fb055fb7ca43e8b854a30786a394a65690a583fbb88c4c850a7fcf78daf75074603484a1c06a742ea4a9d002bf9b63808aeee5759acee12b509649987d7270d3c561273221ebfbba91d5a0c2
    n = 0x26a4159b3400bfed201dba82df5c3f55f75b70984638ccda45d4079e5d3b97c8b78d2bb8fbff84726afe91c6c4112fb96ca0a1716c12a0eae299b835cd4c05623913386752579775193e447b5ebf1b530b78dc7b5bcfedfb337885eae68ea3a4b994ee7ea2a443d8c1daf95bc29d0b37b8037ae7968df83ff7c1a7a9523b6b78042eb44c677662c2

    X = 0x9f910b5ab8f0c95906b67d8867694eae3f74da06d105ba2a05e1a459b29c00a0d109
    Y = 0x13ec07890859f0d2fdd0d79e517fb6f49886c959639a1ef72cc8fd885c000012e
    T = 0xffc00020fffffffc
    NAF_T = NAF(T)

    Zx = ZZ['t']  # Define the ring of polynomials with integer coefficients
    (t,) = Zx._first_ngens(1)
    xt = t ** 4 + 4 * Y * t ** 3 + 8 * Y ** 2 * t ** 2 + 4 * p * Y * t + p ** 2
    assert n == xt(t=1)  # Check is the order is correct

    lx, ly = 41009, 1
    U = [T, (T // 4), lx, ly]

    Fp = GF(p, proof=False)  # Fix the prime field Fp
    Fpx = Fp['x']  # Fix Fpx as the ring of polynomials in x, with coefficients in Fp
    (x,) = Fpx._first_ngens(1)

    # Hyperelliptic curve
    C = HyperellipticCurve(x ** 5 + 3 * x)  # Set the equation of the hyperelliptic curve C
    Jac = C.jacobian()  # J is the Jacobian of the curve C
    # print(Jac)

    P = JC_random_element(C)
    P = h * P  # Force P to have order r
    # Q = JC_random_element(C)
    # h = n // r

    # Constructing the curve and Jacobian over the extension field Fp**8
    # print("Fp12 = Fp[w]/(w^12 + w^6 + 2)")
    Fpw = Fp['w']  # Fpw: polynomial ring in w with coefficients in Fp
    (w,) = Fpw._first_ngens(1)
    Fp8 = Fp.extension(w ** 8 + 3, names=('w',))  # Fp8: degree 8 extension of Fp.
    (w,) = Fp8._first_ngens(1)

    Fp8x = Fp8['x']
    (x,) = Fp8x._first_ngens(1)
    C8 = HyperellipticCurve(Fp8x([0, 3, 0, 0, 0, 1]))
    # print(C8)
    Jac8 = C8.jacobian()
    # print(Jac8)

    Res = (t ** 8 - 1).resultant(xt)  # Compute the resultant  of the polynomials t^8 - 1 and \chi(t)
    h_ = Res // r ** 4  # Compute the cofactor of the Jacobian J8
    n_ = h_ * r  # Compute the order of the Jacobian J8

    Q = JC_random_element(C8)  # Pick random element in Jacobian J8(Fp8)
    Q = h_ * Q  # Force Q to have order r

    print(P)
    print(Q)

    return U


generate_jacobian()
