import pytest

from atom2seq.classes import Atom
from atom2seq.classes import IndexedObject as IO


@pytest.fixture
def atoms() -> tuple[Atom]:
    IO.used_indices = {}
    atom1 = Atom("C", (0, 0, 0))
    atom2 = Atom("H", (0, 1, 0))
    return (atom1, atom2)


def test_symbol(atoms):
    atom1, atom2 = atoms
    assert (atom1.symbol == "C") and (atom2.symbol == "H")


def test_coords(atoms):
    atom1, atom2 = atoms
    assert (atom1.coords == (0, 0, 0)) and (atom2.coords == (0, 1, 0))


def test_repr(atoms):
    atom1, atom2 = atoms
    assert repr(atom1) == "Atom('C', (0, 0, 0), -1, 0)"


def test_get_idx(atoms):
    atom1, atom2 = atoms
    assert (atom1.get_idx() == 0) and (atom2.get_idx() == 1)


def test_set_idx(atoms):
    atom1, atom2 = atoms
    atom1.set_idx(2)
    atom2.set_idx(3)
    assert (atom1.get_idx() == 2) and (atom2.get_idx() == 3)


def test_get_parent(atoms):
    atom1, atom2 = atoms
    assert atom1.get_parent() == atom2.get_parent() == -1


def test_set_parent(atoms):
    atom1, atom2 = atoms
    atom1.set_parent(1)
    atom2.set_parent(2)
    assert (atom1.get_parent() == 1) and (atom2.get_parent() == 2)


def test_eq(atoms):
    atom1, atom2 = atoms
    atom3 = Atom(atom1.symbol, atom1.coords)
    assert (atom1 == atom3) and (atom1 != atom2)


def test_dist(atoms):
    atom1, atom2 = atoms
    assert atom1.dist(atom2) == 1
