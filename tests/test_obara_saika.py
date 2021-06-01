from obarasaika.obara_saika import find_fun_to_lower, find_component_to_lower


def test_find_fun_to_lower():
    assert find_fun_to_lower([1, 0, 0, 0, 0, 0], 2) == 0
    assert find_fun_to_lower([0, 1, 0, 0, 0, 0], 2) == 0
    assert find_fun_to_lower([0, 0, 1, 0, 0, 0], 2) == 0
    assert find_fun_to_lower([0, 0, 0, 1, 0, 0], 2) == 1
    assert find_fun_to_lower([0, 0, 0, 0, 1, 0], 2) == 1
    assert find_fun_to_lower([0, 0, 0, 0, 0, 1], 2) == 1
    assert find_fun_to_lower([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 4) == 1
    assert find_fun_to_lower([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 4) == 3
    assert find_fun_to_lower([1, 0, 0, 0, 0, 0, 0, 0, 1], 3) == 0
    assert find_fun_to_lower([0, 0, 0, 0, 0, 0, 0, 0, 1], 3) == 2
    assert find_fun_to_lower([0, 0, 0, 0, 1, 0, 0, 0, 1], 3) == 1
    assert find_fun_to_lower([0, 0, 0, 0, 2, 0, 0, 0, 1], 3) == 2
    assert find_fun_to_lower([0, 0, 0, 0, 2, 0, 0, 0, 1], 2) == 1


def test_find_component_to_lower():
    assert find_component_to_lower([0, 0, 1]) == 2
    assert find_component_to_lower([0, 1, 1]) == 1
    assert find_component_to_lower([1, 0, 1]) == 0
    assert find_component_to_lower([0, 0, 0]) == -1
