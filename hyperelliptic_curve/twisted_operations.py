from operations import ADD, DBL
from final_exponentiation import final_exponentiation


def Twisted_Ate_k8(D, Dn, Q_vec, Q, F, s, U, K, case: str = 'case1'):
    T, fc = D, 1
#    for i in s:
    for i in range(len(s)-2, -1, -1):
        T, lc = DBL(T, Q_vec, Q, F, case)
        fc = lc * fc**2
        if s[i] == 1:
            T, lc = ADD(D, T, Q_vec, Q, F, case)
            fc = lc * fc
        if s[i] == -1:
            T, lc = ADD(Dn, T, Q_vec, Q, F, case)
            fc = lc * fc
    fc = final_exponentiation(fc, U, K)
    return fc
