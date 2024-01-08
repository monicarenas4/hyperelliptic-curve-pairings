from sage.rings.integer_ring import ZZ

from generate_jacobian_cp8 import generate_jacobian
from generate_jacobian_k16 import generate_jacobian_k16
import datetime
from os.path import exists
from _utils import make_folder, head_operations_file, write_results

TODAY = str(datetime.date.today()).replace('-', '')
make_folder('results')
txt_results = 'results/number_of_operations.txt'
head_operations_file(txt_results) if not exists(txt_results) else None

print('+++++++++++++++++++++\nExample Jacobian CP8\n+++++++++++++++++++++')
generate_jacobian()

print('+++++++++++++++++++++\nExample Jacobian k16\n+++++++++++++++++++++')
u = ZZ(0x100004003)
generate_jacobian_k16(u)
