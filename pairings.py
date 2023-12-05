from sage.all_cmdline import *

from generate_curves import BLS12_curve
from miller_loop import miller_loop_tate_pairing
from final_exponentiation import final_exponentiation_BLS12
from correctness_tests import tate_pairing_bilinearity_check

import datetime
import time
from os.path import exists
from utils_ import make_folder, head_text_file

TODAY = str(datetime.date.today()).replace('-', '')
make_folder('results')

txt_results = 'results' + '/' + TODAY + ' time_pairing.txt'
head_text_file(txt_results) if not exists(txt_results) else None


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
    pairing_value = final_exponentiation_BLS12(miller_function, p, r)
    tf_pairing = round(time.time() - t0_pairing, 6)

    with open(txt_results, 'a') as file:
        file.write(
            f'{tf_miller}' + '\t' +
            f'{tf_pairing}' + '\t' +
            f'{tf_miller + tf_pairing}' +
            '\n')

    # Bilinearity check
    tate_pairing_bilinearity_check(P, Q, pairing_value, p, r) if check == True else None

    return miller_function, pairing_value
