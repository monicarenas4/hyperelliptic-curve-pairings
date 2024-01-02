from sage.all import Integer
from jacobian_operations import ADD, DBL
from _utils import NAF

def miller_function(P, Q, Q_prec, c_vec, F, length_miller, case, twist = None):
    length_miller = Integer(length_miller).digits(2)
    T, fc = P, 1

    for i in range(len(length_miller) - 2, -1, -1):
        T, lc = DBL(T, Q_prec, Q, F, c_vec, case, twist)
        fc = lc * fc ** 2
        if length_miller[i] == 1:
            T, lc = ADD(P, T, Q_prec, Q, F, c_vec, case, twist)
            fc = lc * fc

    return fc

def miller_function_naf(P, P_neg, Q, Q_prec, c_vec, F, length_miller, case, twist = None):
    length_miller_naf = NAF(length_miller)
    T, fc = P, 1

    for i in range(len(length_miller_naf) - 2, -1, -1):
        T, lc = DBL(T, Q_prec, Q, F, c_vec, case, twist)
        fc = lc * fc ** 2
        if length_miller_naf[i] == 1:
            T, lc = ADD(P, T, Q_prec, Q, F, c_vec, case, twist)
            fc = lc * fc
        elif length_miller_naf[i] == -1:
            T, lc = ADD(P_neg, T, Q_prec, Q, F, c_vec, case, twist)
            fc = lc * fc
            
    return fc
