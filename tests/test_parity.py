import itertools

import pytest

from rxnchirality.parity import descriptor, permutation_parity, same_handedness


def test_identity_and_single_swap():
    assert permutation_parity("abcd", "abcd") == 1
    assert permutation_parity("abcd", "bacd") == -1
    assert permutation_parity("abcd", "bcad") == 1  # 3-cycle
    assert permutation_parity("abcd", "dabc") == -1  # 4-cycle


def test_parity_is_multiplicative():
    base = "abcd"
    for p in itertools.permutations(base):
        for q in itertools.permutations(base):
            assert permutation_parity(base, q) == permutation_parity(base, p) * permutation_parity(
                p, q
            )


def test_not_a_permutation():
    with pytest.raises(ValueError):
        permutation_parity("abcd", "abce")
    with pytest.raises(ValueError):
        permutation_parity("aacd", "acad")


def test_same_handedness_any_listing():
    tokens, sign = (1, 2, 3, 4), 1
    for p in itertools.permutations(tokens):
        par = permutation_parity(tokens, p)
        assert same_handedness(tokens, sign, p, sign * par)
        assert not same_handedness(tokens, sign, p, -sign * par)


def test_descriptor():
    rank = {"a": 3, "b": 2, "c": 1, "d": 0}
    # Lowest first, then highest->lowest, CCW seen from the lowest = R.
    assert descriptor(("d", "a", "b", "c"), 1, rank) == "R"
    assert descriptor(("d", "a", "b", "c"), -1, rank) == "S"
    assert descriptor(("a", "d", "b", "c"), 1, rank) == "S"
    assert descriptor(("a", "b", "c", "d"), 1, {"a": 1, "b": 1, "c": 2, "d": 0}) is None
