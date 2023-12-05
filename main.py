from sage.all_cmdline import *
from sage.rings.integer_ring import ZZ

from pairings import tate_pairing

# Tate Pairing
u = -ZZ(2 ** 63 + 2 ** 62 + 2 ** 60 + 2 ** 57 + 2 ** 48 + 2 ** 16)
tate_pairing(u, check=True)
