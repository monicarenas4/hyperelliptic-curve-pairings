def operations_main(file_name: str, type: str):
    with open(file_name, 'a') as file:
        file.write(type + '\n')
    return


def write_number_operations_head(file_name: str):
    with open(file_name, 'a') as file:
        file.write('function' + '\t' + 'embedding_degree' + '\t'
                   + 'case' + '\t' + 'twist' + '\t'
                   + 'NAF_rep' + '\t'
                   + 'mult_pre' + '\t' + 'sq_pre' + '\t'
                   + 'number_of_DBL' + '\t'
                   + 'mult_DBL' + '\t' + 'sq_DBL' + '\t'
                   + 'mult_miller_DBL' + '\t' + 'sq_miller_DBL' + '\t'
                   + 'number_of_ADD' + '\t'
                   + 'mult_ADD' + '\t' + 'sq_ADD' + '\t'
                   + 'mult_miller_ADD' + '\t' + 'sq_miller_ADD' + '\t'
                   + 'exp_u' + '\t' + 'exp_u0' + '\t'
                   + 'exp_lx' + '\t' + 'exp_um' + '\t' + 'exp_up' + '\t'
                   + 'mult_FE' + '\t' + 'sq_FE' + '\t' + 'inv_FE' + '\t'
                   + 'frobenius' + '\t'
                   + 'total'
                   + '\n')
    return


def write_number_operations(file_name: str = None, function: str = None,
                            embedding_degree: int = None,
                            case: str = None, twist=None,
                            NAF_rep: bool = False,
                            mult_pre: int = 0, sq_pre: int = 0,
                            number_of_DBL: int = 0,
                            mult_DBL: int = 0, sq_DBL: int = 0,
                            mult_miller_DBL: int = 0, sq_miller_DBL: int = 0,
                            number_of_ADD: int = 0,
                            mult_ADD: int = 0, sq_ADD: int = 0,
                            mult_miller_ADD: int = 0, sq_miller_ADD: int = 0,
                            exp_u: int = 0, exp_u0: int = 0, exp_lx: int = 0,
                            exp_um: int = 0, exp_up: int = 0,
                            mult_FE: int = 0, sq_FE: int = 0, inv_FE: int = 0,
                            frobenius: int = 0,
                            total: int = 0):
    with open(file_name, 'a') as file:
        file.write(function + '\t'
                   + str(embedding_degree) + '\t'
                   + case + '\t'
                   + f'{twist}' + '\t'
                   + f'{NAF_rep}' + '\t'
                   + f'{mult_pre}' + '\t' + f'{sq_pre}' + '\t'
                   + f'{number_of_DBL}' + '\t'
                   + f'{mult_DBL}' + '\t' + f'{sq_DBL}' + '\t'
                   + f'{mult_miller_DBL}' + '\t' + f'{sq_miller_DBL}' + '\t'
                   + f'{number_of_ADD}' + '\t'
                   + f'{mult_ADD}' + '\t' + f'{sq_ADD}' + '\t'
                   + f'{mult_miller_ADD}' + '\t' + f'{sq_miller_ADD}' + '\t'
                   + f'{exp_u}' + '\t'
                   + f'{exp_u0}' + '\t'
                   + f'{exp_lx}' + '\t'
                   + f'{exp_um}' + '\t'
                   + f'{exp_up}' + '\t'
                   + f'{mult_FE}' + '\t' + f'{sq_FE}' + '\t' + f'{inv_FE}' + '\t'
                   + f'{frobenius}' + '\t'
                   + f'{total}'
                   + '\n')
    return
