from operations import ADD, DBL
from final_exponentiation import final_exponentiation


def Twisted_Ate_k8(D, Dn, E, F, s, U, K, pow):
    T, fc = D, 1

    for i in range(len(s)-2, -1, -1):
        T, lc = DBL(T, E, F)
        fc = lc * fc**2
        if s[i] == 1:
            T, lc = ADD(D, T, E, F)
            fc = lc * fc
        if s[i] == -1:
            T, lc = ADD(Dn, T, E, F)
            fc = lc * fc
    fc = final_exponentiation(fc, U, K)
    
    return fc
