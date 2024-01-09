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


def operations_main(file_name: str, type: str):
    with open(file_name, 'a') as file:
        file.write(type
                   + '\n')
    return


def operations_bilinearity_check_head(file_name):
    with open(file_name, 'a') as file:
        file.write('function' + '\t'
                   + 'case' + '\t'
                   + 'NAF_rep' + '\t'
                   + 'mult_pre' + '\t'
                   + 'sq_pre' + '\t'
                   + '\n')
    return


def operations_bilinearity_check(file_name, function, case, NAF_rep, mult_pre, sq_pre):
    with open(file_name, 'a') as file:
        file.write(function + '\t'
                   + case + '\t'
                   + f'{NAF_rep}' + '\t'
                   + f'{mult_pre}' + '\t'
                   + f'{sq_pre}' + '\t'
                   + '\n')
    return


def operations_miller_loop_head(file_name: str):
    with open(file_name, 'a') as file:
        file.write('function' + '\t'
                   + 'case' + '\t'
                   + 'twist' + '\t'
                   + 'NAF_rep' + '\t'
                   + 'size_DBL' + '\t'
                   + 'size_ADD' + '\t'
                   + 'mult_line' + '\t'
                   + 'sq_line' + '\t'
                   + 'mult_DBL' + '\t'
                   + 'sq_DBL' + '\t'
                   + 'mult_ADD' + '\t'
                   + 'sq_ADD' + '\t'
                   + 'mult_miller' + '\t'
                   + 'sq_miller' + '\t'
                   + '\n')
    return


def operations_miller_loop(file_name: str, function: str, case: str, twist, NAF_rep: bool,
                           size_DBL, size_ADD,
                           mult_line: int, sq_line: int,
                           mult_DBL: int, sq_DBL: int,
                           mult_ADD: int, sq_ADD: int,
                           mult_miller: int, sq_miller: int):
    with open(file_name, 'a') as file:
        file.write(function + '\t'
                   + case + '\t'
                   + f'{twist}' + '\t'
                   + f'{NAF_rep}' + '\t'
                   + f'{size_DBL}' + '\t'
                   + f'{size_ADD}' + '\t'
                   + f'{mult_line}' + '\t'
                   + f'{sq_line}' + '\t'
                   + f'{mult_DBL}' + '\t'
                   + f'{sq_DBL}' + '\t'
                   + f'{mult_ADD}' + '\t'
                   + f'{sq_ADD}' + '\t'
                   + f'{mult_miller}' + '\t'
                   + f'{sq_miller}' + '\t'
                   + '\n')
    return
