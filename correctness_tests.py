from sage.all_cmdline import *

from pairings import miller_loop_tate_pairing, miller_loop_opt_tate_pairing, miller_loop_opt_ate_pairing
from pairings import final_exp_BLS12, final_exp_BLS12_optimized
from random import randint


def pairing_bilinearity_check(P, Q, pairing_value, p: int, r: int, type: str = 'tate', u=None):
    """
    :param P: P-point in (x,y,z)-coordinates
    :param Q: Q-point in (x,y,z)-coordinates
    :param pairing_value:
    :param p:
    :param r:
    :param type: str
    :param u: int
    :return:
    """
    global pairing_value_prime
    a, b = randint(0, r - 1), randint(0, r - 1)
    pairing_value = pairing_value ** (a * b)
    # print('Pairing value 1:', pairing_value)

    P_prime, Q_prime = a * P, b * Q  # P' and Q' pairings
    if 'tate' == type:
        miller_function_prime = miller_loop_tate_pairing(P_prime, Q_prime, r)
        pairing_value_prime = final_exp_BLS12(miller_function_prime, p, r)
    elif 'optimized_tate' == type:
        miller_function_prime = miller_loop_opt_tate_pairing(P_prime, Q_prime, r)
        pairing_value_prime = final_exp_BLS12_optimized(miller_function_prime, u, p)
    elif 'optimal_ate' == type:
        miller_function_prime = miller_loop_opt_ate_pairing(Q_prime, P_prime, u)
        pairing_value_prime = final_exp_BLS12_optimized(miller_function_prime, u, p)

    # print('Pairing value 2:', pairing_value_prime)
    print("Bilinearity Check:", (pairing_value == pairing_value_prime))

    return None
