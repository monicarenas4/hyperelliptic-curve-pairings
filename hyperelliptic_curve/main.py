from sage.rings.integer_ring import ZZ
import datetime

from generate_jacobian_cp8 import generate_jacobian
from generate_jacobian_k16 import generate_jacobian_k16
from generate_jacobian_k24 import generate_jacobian_k24
from _utils import make_folder
from write_number_operations import operations_main, write_number_operations_head

TODAY = str(datetime.date.today()).replace('-', '')
make_folder('results')
file_name = 'results/number_of_operations.txt'

write_number_operations_head(file_name)

print('+++++++++++++++++++++\nExample Jacobian CP8\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian CP8')
generate_jacobian()

print('+++++++++++++++++++++\nExample Jacobian Kawazoe-Takahashi k = 16\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian k16')
u = ZZ(0x100004003)
assert u == 2 ** 32 + 2 ** 14 + 2 ** 1 + 2 ** 0
generate_jacobian_k16(u, family='k16')

print('+++++++++++++++++++++\nExample Jacobian from new family for k = 16\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian new family k16')
u = ZZ(0xeffeffff7fff)
generate_jacobian_k16(u, family='new_k16')

print('+++++++++++++++++++++\nExample Jacobian from new family for k = 24\nFirst seed\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian family k24')
u0 = ZZ(49)
u1 = ZZ(2**48 + 2**32 + 2**0)
u2 = ZZ(2**47 + 2**34 + 2**4 + 2**0)
u3 = ZZ(2**10 - 2**7 + 2**5 + 2**3 + 1)
u4 = ZZ(2**48 + 2**44 + 2**6 + 2**4 + 2**3 + 1)
generate_jacobian_k24(u1, seed="u1", family='k24')

print('+++++++++++++++++++++\nExample Jacobian from new family for k = 24\nSecond seed\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian family k24')
generate_jacobian_k24(u4, seed="u4", family='k24')
