from sage.all_cmdline import *

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF

R = QQ['x']
(x,) = R._first_ngens(1)

rx = x ** 4 - x ** 2 + 1
# Trace polynomial
tx = x + 1
px = (x ** 6 - 2 * x ** 5 + 2 * x ** 3 + x + 1) / 3

# seed
u = -ZZ(2 ** 63 + 2 ** 62 + 2 ** 60 + 2 ** 57 + 2 ** 48 + 2 ** 16)

r = ZZ(rx(u))
t = ZZ(tx(u))
p = ZZ(px(u))
n = p - t + 1  # Order of the curve. #E = p-t+1 = hr
h = n / r  # Co-factor

Fp = GF(p, proof=False)
E = EllipticCurve([Fp(0), Fp(-3)])
P = E.random_element()
P = h * P

print("Random Point:", P)
print((P[1] ** 2 - P[0] ** 3 + 3) % p)
print(r * P)
# print('The order of the curve is:', n)
# print(n%r)
# print((p**12-1)%r)
# print(h)

# print(r, is_prime(r))
# print(t, is_prime(t))
# print (p, is_prime(p))

Fpw = Fp['w']; (w,) = Fpw._first_ngens(1)

print("Fp12 = Fp[w]/(w^12 + w^6 + 2)")
Fp12 = Fp.extension(w**12 + w**6 + 2, names=('w',)); (w,) = Fp12._first_ngens(1)
print(Fp12)
E12 = EllipticCurve([Fp12(0), Fp12(-3)])
Q = E12.random_element()
print(E12.order()%r)





