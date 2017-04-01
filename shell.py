def get_ijk_list(m):
    """Form all possible (i, j, k) exponents up to maximum total angular
    momentum m.
    """

    l = []
    for a in range(1, m + 2):
        for b in range(1, a + 1):
            i = m + 1 - a
            j = a - b
            k = b - 1
            l.append([i, j, k])
    return l


def get_shell4(a, b, c, d):
    """Form all possible angular momentum combinations up to lmax = a, b,
    c, d for center 1, 2, 3, and 4, respectively.
    """

    components = []
    for p in get_ijk_list(a):
        for q in get_ijk_list(b):
            for r in get_ijk_list(c):
                for s in get_ijk_list(d):
                    components.append(p + q + r + s)
    return components


def get_shell2(a, b):
    """Form all possible angular momentum combinations up to lmax = a on
    center 1 and lmax = b on center 2.
    """

    components = []
    for p in get_ijk_list(a):
        for q in get_ijk_list(b):
            components.append(p + q)
    return components
