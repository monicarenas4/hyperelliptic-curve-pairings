def final_exponentiation_BLS12(miller_fn, p, r, k: int = 12) -> float:
    """
    :param miller_fn: miller function
    :param p: prime number
    :param r: prime number
    :param k: embedding factor
    :return: exponentiation
    """
    exponentiation = miller_fn ** ((p ** k - 1) / r)

    return exponentiation
