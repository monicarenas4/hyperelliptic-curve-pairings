from sage.rings.integer_ring import ZZ

from generate_jacobian_cp8 import generate_jacobian
from generate_jacobian_k16 import generate_jacobian_k16
from generate_jacobian_new_k16 import generate_jacobian_new_k16
import datetime
from os.path import exists
from _utils import make_folder, head_operations_file, write_results

TODAY = str(datetime.date.today()).replace('-', '')
make_folder('results')
txt_results = 'results/number_of_operations.txt'
head_operations_file(txt_results) if not exists(txt_results) else None

print('Example Jacobian CP8')
generate_jacobian()

print('Example Jacobian Kawazoe-Takahashi k = 16')
u = ZZ(0x100004003)
generate_jacobian_k16(u, family='k16')

print('Example Jacobian from new family for k = 16')
u = ZZ(0xeffeffff7fff)
generate_jacobian_new_k16(u, family='afk16')
