import pytest

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


@pytest.fixture
def hydroxyl():
    Atom.idxs = {}
    atoms = set()
    atoms.add(Atom("O", (0, 0, 0)))
    atoms.add(Atom("H", (0, 0, 1)))
    bonds = ConnectivityTable({(0, 1)})
    return Group(atoms, bonds)


def test_repr(hydroxyl):
    assert (
        repr(hydroxyl) == "Group([Atom('O', (0, 0, 0), 0), Atom('H', "
        "(0, 0, 1), 1)], ConnectivityTable({(0, 1)}))"
    )


def test_lt_len_atoms(hydroxyl):
    assert hydroxyl > Group(set(), ConnectivityTable(set()))


def test_lt_len_bonds(hydroxyl):
    assert hydroxyl > Group(
        {Atom("O", (0, 0, 0), 0), Atom("H", (0, 0, 1), 1)},
        ConnectivityTable(set()),  # noqa
    )


def test_lt_by_atom(hydroxyl):
    assert hydroxyl > Group(
        {Atom("O", (0, 0, -1), 0), Atom("H", (0, 0, 1), 1)},
        ConnectivityTable({(0, 1)}),  # noqa
    )


def test_eq(hydroxyl):
    assert hydroxyl == Group(
        {Atom("O", (0, 0, 0), 2), Atom("H", (0, 0, 1), 3)},
        ConnectivityTable([(0, 1)]),  # noqa
    )


def test_dist(hydroxyl):
    assert hydroxyl.dist(0, 1) == 1


def test_get_bonds(hydroxyl):
    assert hydroxyl.get_bonds() == ConnectivityTable([(0, 1)])


def test_get_atoms(hydroxyl):
    assert hydroxyl.get_atoms() == {Atom("O", (0, 0, 0)), Atom("H", (0, 0, 1))}


def test_atom_list(hydroxyl):
    assert hydroxyl.atom_list() == [Atom("O", (0, 0, 0)), Atom("H", (0, 0, 1))]


def test_get_idx(hydroxyl):
    assert hydroxyl.get_idx() == -1


def test_set_idx(hydroxyl):
    hydroxyl.set_idx(3)
    assert hydroxyl.get_idx() == 3


def test_atom_idxs(hydroxyl):
    assert (hydroxyl.atom_list()[0].get_idx() == 0) and (
        hydroxyl.atom_list()[1].get_idx() == 1
    )
