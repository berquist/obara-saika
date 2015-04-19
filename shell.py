def get_ijk_list(m):
    """
    Form all possible (i, j, k) exponents up to maximum total angular momentum m.
    """
    l = []
    for a in range(1, m + 2):
        for b in range(1, a + 1):
            i = m + 1 - a
            j = a - b
            k = b - 1
            l.append([i, j, k])
    return l


def test_get_ijk_list():
    """
    Tests for get_ijk_list() up to m = 2.
    """

    l0 = get_ijk_list(0)
    assert l0 == [
        [0, 0, 0]
    ]

    l1 = get_ijk_list(1)
    assert l1 == [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]

    l2 = get_ijk_list(2)
    assert l2 == [
        [2, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [0, 2, 0],
        [0, 1, 1],
        [0, 0, 2]
    ]

#-------------------------------------------------------------------------------

def get_shell4(a, b, c, d):
    """
    Form all possible angular momentum combinations up to lmax = a, b, c, d
    for center 1, 2, 3, and 4, respectively.
    """
    components = []
    for p in get_ijk_list(a):
        for q in get_ijk_list(b):
            for r in get_ijk_list(c):
                for s in get_ijk_list(d):
                    components.append(p + q + r + s)
    return components

#-------------------------------------------------------------------------------

def get_shell2(a, b):
    """
    Form all possible angular momentum combinations up to lmax = a on center 1
    and lmax = b on center 2.
    """
    components = []
    for p in get_ijk_list(a):
        for q in get_ijk_list(b):
            components.append(p + q)
    return components


def test_get_shell2():
    """
    Some low angular momentum tests for get_shell2().
    """

    a0b0 = get_shell2(0, 0)
    assert a0b0 == [
        [0, 0, 0, 0, 0, 0]
    ]

    a0b1 = get_shell2(0, 1)
    a1b0 = get_shell2(1, 0)
    assert a0b1 == [
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1]
    ]
    assert a1b0 == [
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0]
    ]

    a0b2 = get_shell2(0, 2)
    assert a0b2 == [
        [0, 0, 0, 2, 0, 0],
        [0, 0, 0, 1, 1, 0],
        [0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 2, 0],
        [0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 2]
    ]
