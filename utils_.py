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
        f.write('miller_loop' + '\t'
                + 'final_exp' + '\t'
                + 'total'
                + '\n')
    return
