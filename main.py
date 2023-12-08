from sage.all_cmdline import *
from sage.rings.integer_ring import ZZ

from pairings import tate_pairing, optimized_tate_pairing
import datetime
from os.path import exists
from utils_ import make_folder, head_text_file, write_results

TODAY = str(datetime.date.today()).replace('-', '')
make_folder('results')
txt_results = 'results/time_pairing.txt'
head_text_file(txt_results) if not exists(txt_results) else None

# Tate Pairing
u = -ZZ(2 ** 63 + 2 ** 62 + 2 ** 60 + 2 ** 57 + 2 ** 48 + 2 ** 16)
_, _, time_miller, time_pairing = tate_pairing(u, check=True)
write_results(TODAY, txt_results, 'tate_pairing', time_miller, time_pairing)

# Optimized Tate Pairing
_, _, tf_miller, tf_total = optimized_tate_pairing(u, check=True)
write_results(TODAY, txt_results, 'opt_tate_pairing', tf_miller, tf_total)
