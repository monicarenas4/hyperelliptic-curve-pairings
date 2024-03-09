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

print('+++++++++++++++++++++\nExample Jacobian Cocks-Pinch k = 8: CP8-544 Curve\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian CP8')
generate_jacobian()

print('\n')
print('+++++++++++++++++++++\nExample Jacobian Kawazoe-Takahashi k = 16: KT16-447 Curve\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian KT16')
# Seed u to generate the KT16-447 curve
u = ZZ(0x100004003)
assert u == 2 ** 32 + 2 ** 14 + 2 ** 1 + 2 ** 0
generate_jacobian_k16(u, family='KT16')

print('\n')
print('+++++++++++++++++++++\nExample Jacobian from new family for k = 16: New16-767\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian new family New16')
u = ZZ(0xeffeffff7fff)
assert u == 2 ** 48 - 2 ** 44 - 2 ** 32 - 2 ** 15 - 1
generate_jacobian_k16(u, family='New16')

print('\n')
print('+++++++++++++++++++++\nExample Jacobian from new family for k = 24\nFirst seed: Curve New24-575\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian family New24')
u1 = ZZ(0x1000100000001)
assert u1 == 2 ** 48 + 2 ** 32 + 2 ** 0
generate_jacobian_k24(u1, seed="u1", family='New24')

print('\n')
print('+++++++++++++++++++++\nExample Jacobian from new family for k = 24\nSecond seed: Curve New24-576\n+++++++++++++++++++++')
operations_main(file_name, 'Jacobian family New24')
u2 = ZZ(2**48 + 2**44 + 2**6 + 2**4 + 2**3 + 1)
generate_jacobian_k24(u2, seed="u2", family='New24')
