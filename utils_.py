import os


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
