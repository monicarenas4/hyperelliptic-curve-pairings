from sage.all_cmdline import *

from elliptic_curve.generate_curves import BLS12_curve
from elliptic_curve.miller_loop import miller_loop_tate_pairing, miller_loop_opt_tate_pairing, \
    miller_loop_opt_ate_pairing
from elliptic_curve.final_exponentiation import final_exp_BLS12, final_exp_BLS12_optimized
from elliptic_curve.correctness_tests import pairing_bilinearity_check
import time


def tate_pairing(u, check=False):
    """
    :param u: curve seed
    :return: m, pairing value
    """
    P, Q, p, r = BLS12_curve(u)
    t0_miller = time.time()
    miller_function = miller_loop_tate_pairing(P, Q, r)
    tf_miller = round(time.time() - t0_miller, 6)

    t0_pairing = time.time()
    pairing_value = final_exp_BLS12(miller_function, p, r)
    tf_pairing = round(time.time() - t0_pairing, 6)

    # Bilinearity check
    pairing_bilinearity_check(P, Q, pairing_value, p, r, type='tate') if check == True else None

    return miller_function, pairing_value, tf_miller, tf_pairing


def optimized_tate_pairing(u, check=False):
    """
    :param u: curve seed
    :return: m, pairing value
    """
    P, Q, p, r = BLS12_curve(u)
    t0_miller = time.time()
    miller_function = miller_loop_opt_tate_pairing(P, Q, r)
    tf_miller = round(time.time() - t0_miller, 6)

    t0_total = time.time()
    pairing_value = final_exp_BLS12_optimized(miller_function, u, p)
    tf_total = round(time.time() - t0_total, 6)

    # Bilinearity check
    pairing_bilinearity_check(P, Q, pairing_value, p, r, type='optimized_tate', u=u) if check == True else None

    return miller_function, pairing_value, tf_miller, tf_total


def optimal_ate_pairing(u, check=False):
    """
    :param u: curve seed
    :return: m, pairing value
    """
    P, Q, p, r = BLS12_curve(u, type='ate')
    t0_miller = time.time()
    miller_function = miller_loop_opt_ate_pairing(Q, P, u)
    tf_miller = round(time.time() - t0_miller, 6)

    t0_total = time.time()
    pairing_value = final_exp_BLS12_optimized(miller_function, u, p)
    tf_total = round(time.time() - t0_total, 6)

    # Bilinearity check
    pairing_bilinearity_check(P, Q, pairing_value, p, r, type='optimal_ate', u=u) if check == True else None

    return miller_function, pairing_value, tf_miller, tf_total
