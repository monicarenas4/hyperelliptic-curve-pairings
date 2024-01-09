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

print('Example Jacobian CP8')
generate_jacobian()

print('Example Jacobian k16')
u = ZZ(0x100004003)
generate_jacobian_k16(u)
