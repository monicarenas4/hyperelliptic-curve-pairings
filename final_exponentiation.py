from sage.all import Integer


def final_exp_BLS12(miller_fn, p, r, k: int = 12):
    """
    :param miller_fn: miller function
    :param p: prime number
    :param r: prime number
    :param k: embedding factor
    :return: exponentiation
    """
    exponentiation = miller_fn ** ((p ** k - 1) / r)

    return exponentiation


def final_exp_easy_BLS12(miller_fn, p: int):
    """
    :param miller_fn: miller function
    :param p: prime number
    :return: exponentiation
    """
    # t0_easy = time.time()
    inverse = miller_fn ** (-1)
    f = miller_fn ** (p ** 6)
    f = f * inverse
    f2 = f ** (p ** 2)
    easy_k12 = f2 * f
    # tf_easy = round(time.time() - t0_easy, 6)

    return easy_k12


def final_exp_hard_BLS12(input_fn, u, p):
    u0 = -u
    # hard_BLS12 = miller_fn ** ((p ** 4 - p ** 2 + 1) / r)
    f3 = input_fn ** (u0 + 1)
    f3 = f3 ** (u0 + 1)
    # res = f3 ** p  # Frobenius prop.
    res = f3.frobenius(1)
    f2 = f3 ** u0
    # f2 = 1 / f2 # Frobenius prop.
    f2 = f2.frobenius(6)
    res = res * f2
    f1 = f2 ** u0
    f1 = f1 * f3
    f1 = 1 / f1
    res = res ** p
    res = res * f1
    f0 = f1 ** u0
    # f0 = 1 / f0 # Frobenius prop.
    f0 = f0.frobenius(6)
    f0 = f0 * input_fn ** 2 * input_fn
    res = res ** p
    hard_BLS12 = res * f0

    return hard_BLS12


def final_exp_BLS12_optimized(miller_fn, u, p: int):
    exp_easy = final_exp_easy_BLS12(miller_fn, p)
    pairing_value = final_exp_hard_BLS12(exp_easy, u, p)

    return pairing_value
