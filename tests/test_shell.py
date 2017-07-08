from obarasaika.shell import get_ijk_list, get_shell2

def test_get_ijk_list():
    """Tests for get_ijk_list() up to m = 2."""

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
