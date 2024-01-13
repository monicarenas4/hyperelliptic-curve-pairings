from sage.rings.integer_ring import ZZ

from generate_jacobian_cp8 import generate_jacobian
from generate_jacobian_k16 import generate_jacobian_k16
import datetime
from _utils import make_folder
from write_number_operations import operations_main, write_number_operations_head

TODAY = str(datetime.date.today()).replace('-', '')
make_folder('results')
file_name = 'results/number_of_operations.txt'

print('+++++++++++++++++++++\nExample Jacobian CP8\n+++++++++++++++++++++')
write_number_operations_head(file_name)
operations_main(file_name, 'Jacobian CP8')
generate_jacobian()

print('+++++++++++++++++++++\nExample Jacobian Kawazoe-Takahashi k = 16\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian k16')
u = ZZ(0x100004003)
generate_jacobian_k16(u, family='k16')

print('+++++++++++++++++++++\nExample Jacobian from new family for k = 16\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian new family k16')
u = ZZ(0xeffeffff7fff)
generate_jacobian_k16(u, family='new_k16')
