from sage.all_cmdline import *

from pairings import miller_loop_tate_pairing, final_exp_BLS12, final_exp_BLS12_optimized
from random import randint


def tate_pairing_bilinearity_check(P, Q, pairing_value, p: int, r: int, type: str = 'tate', seed=None):
    a, b = randint(0, r - 1), randint(0, r - 1)
    pairing_value = pairing_value ** (a * b)

    # pairing of P' and Q'
    P_prime, Q_prime = a * P, b * Q
    miller_function_prime = miller_loop_tate_pairing(P_prime, Q_prime, r)
    if not 'optimized' in type:
        pairing_value_prime = final_exp_BLS12(miller_function_prime, p, r)
    else:
        pairing_value_prime = final_exp_BLS12_optimized(miller_function_prime, seed, p)

    print("Bilinearity Check:", (pairing_value == pairing_value_prime))

    return None
