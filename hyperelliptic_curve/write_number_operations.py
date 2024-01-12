def operations_main(file_name: str, type: str):
    with open(file_name, 'a') as file:
        file.write(type + '\n')
    return


def write_number_operations_head(file_name: str):
    with open(file_name, 'a') as file:
        file.write('function' + '\t'
                   + 'case' + '\t' + 'twist' + '\t'
                   + 'NAF_rep' + '\t'
                   + 'mult_pre' + '\t' + 'sq_pre' + '\t'
                   + 'size_DBL' + '\t' + 'size_ADD' + '\t'
                   + 'mult_line' + '\t' + 'sq_line' + '\t'
                   + 'mult_DBL' + '\t' + 'sq_DBL' + '\t'
                   + 'mult_ADD' + '\t' + 'sq_ADD' + '\t'
                   + 'mult_miller' + '\t' + 'sq_miller' + '\t'
                   + 'exp_u' + '\t' + 'exp_u0' + '\t'
                   + 'exp_lx' + '\t' + 'exp_um' + '\t'
                   + 'mult_FE' + '\t' + 'sq_FE' + '\t' + 'inv_FE' + '\t'
                   + 'frobenius' + '\t' + 'frob_power'  # 25
                   + '\n')
    return


def write_number_operations(file_name: str = None, function: str = None, case: str = None, twist=None,
                            NAF_rep: bool = None,
                            mult_pre: int = None, sq_pre: int = None,
                            size_DBL=None, size_ADD=None,
                            mult_line: int = None, sq_line: int = None,
                            mult_DBL: int = None, sq_DBL: int = None,
                            mult_ADD: int = None, sq_ADD: int = None,
                            mult_miller: int = None, sq_miller: int = None,
                            exp_u: int = None, exp_u0: int = None, exp_lx: int = None, exp_um: int = None,
                            mult_FE: int = None, sq_FE: int = None, inv_FE: int = None,
                            frobenius: int = None, frob_power: int = None,
                            factor=1):
    with open(file_name, 'a') as file:
        file.write(function + '\t'
                   + case + '\t'
                   + f'{twist}' + '\t'
                   + f'{NAF_rep}' + '\t'
                   + f'{mult_pre}' + '\t' + f'{sq_pre}' + '\t'
                   + f'{size_DBL if size_DBL is None else size_DBL * factor}' + '\t'
                   + f'{size_ADD if size_ADD is None else size_ADD * factor}' + '\t'
                   + f'{mult_line}' + '\t' + f'{sq_line}' + '\t'
                   + f'{mult_DBL}' + '\t' + f'{sq_DBL}' + '\t'
                   + f'{mult_ADD}' + '\t' + f'{sq_ADD}' + '\t'
                   + f'{mult_miller}' + '\t' + f'{sq_miller}' + '\t'
                   + f'{exp_u}' + '\t' + f'{exp_u0}' + '\t'
                   + f'{exp_lx}' + '\t' + f'{exp_um}' + '\t'
                   + f'{mult_FE}' + '\t' + f'{sq_FE}' + '\t' + f'{inv_FE}'
                   + '\t' + f'{frobenius}' + '\t' + f'{frob_power}'
                   + '\n')
    return
