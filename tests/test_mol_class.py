import pytest

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.mol_class import Mol


@pytest.fixture
def water():
    Atom.idxs = {}
    return Mol(
        {Atom("O", (0, 0, 0)), Atom("H", (0, 1, 0)), Atom("H", (1, 0, 0))},
        ConnectivityTable({(0, 1), (0, 2)}),
    )


def test_repr(water):
    assert (
        repr(water) == "Mol([Atom('O', (0, 0, 0), 0), Atom('H', (0, 1, 0), "
        "1), Atom('H', (1, 0, 0), 2)], ConnectivityTable([(0, 1), (0, 2)]))"
    )


def test_eq(water):
    assert water == Mol(
        {
            Atom("O", (0, 0, 0), 3),
            Atom("H", (0, 1, 0), 4),
            Atom("H", (1, 0, 0), 5),
        },  # noqa
        ConnectivityTable({(0, 1), (0, 2)}),
    )


def test_dist(water):
    assert water.dist(0, 1) == 1


def test_get_bonds(water):
    assert water.get_bonds() == ConnectivityTable({(0, 1), (0, 2)})


def test_get_atoms(water):
    assert water.get_atoms() == {
        Atom("O", (0, 0, 0)),
        Atom("H", (0, 1, 0)),
        Atom("H", (1, 0, 0)),
    }


def test_group_atoms(water):
    hydroxyl = water.group_atoms([0, 1])
    assert hydroxyl == Group(
        {Atom("O", (0, 0, 0), 3), Atom("H", (0, 1, 0), 4)},
        ConnectivityTable({(0, 1)}),  # noqa
    )  # noqa
