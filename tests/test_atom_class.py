import pytest

from atom2seq.atom_class import Atom


@pytest.fixture
def hydrogen():
    return Atom("H", (0, 0, 0))


def test_eq(hydrogen):
    assert hydrogen == Atom("H", (0, 0, 0))


def test_get_idx(hydrogen):
    assert hydrogen.get_idx() == -1


def test_set_idx(hydrogen):
    hydrogen.set_idx(4)
    assert hydrogen.get_idx() == 4
