import os
import time
from sage.all import Integer


def make_folder(folder_name: str):
    """
    :param folder_name: folder name
    :return:
    """
    os.makedirs(folder_name) if not os.path.exists(folder_name) else None
    return


def head_text_file(file_name: str):
    """
    :param file_name: file name
    :return: txt file
    """
    with open(file_name, 'w') as f:
        f.write('date' + '\t'
                + 'pairing_type' + '\t'
                + 'miller_loop' + '\t'
                + 'final_exp' + '\t'
                + 'total_pairing'
                + '\n')
    return


def write_results(date, file_name: str, pairing_name: str, tf_miller: float, tf_pairing: float):
    with open(file_name, 'a') as file:
        file.write(date + '\t'
                   + pairing_name + '\t'
                   + f'{tf_miller}' + '\t'
                   + f'{tf_pairing}' + '\t'
                   + f'{tf_miller + tf_pairing}'
                   + '\n')
    return


def NAF(x: int):
    """
    :param x: int
    :return: list
    """
    t0 = time.time()

    naf_x = []
    xx = Integer(x)
    assert x >= 0
    while xx > 0:
        rr = xx % 4
        if rr == 3:
            rr = -1
        else:
            rr = rr % 2
        naf_x.append(rr)
        xx -= rr
        xx, rr = xx.quo_rem(2)
        assert rr == 0
    assert x == sum([r * 2 ** i for i, r in enumerate(naf_x)])

#    print(time.time() - t0)
    # naf_x1 = naf_x.reverse()

    return naf_x


def hamming_weight(bit_x: list) -> int:
    """
    :param bit_x: binary representartion of a positive integer
    :return: Hamming weight
    """
    count = 0
    for i in range(len(bit_x)):
        if bit_x[i] == 1:
            count = count + 1

    return count


def NAf_hamming_weight(naf_x: list) -> int:
    """
    :param naf_x: NAF representation of a positive integer
    :return: NAF Hamming weight
    """
    count = 0
    for i in range(0, len(naf_x)):
        if (naf_x[i] == 1) or (naf_x[i] == -1):
            count = count + 1

    return count
