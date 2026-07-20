import pytest

from atom2seq.atom_class import Atom
from atom2seq.connect_groups import connect_groups
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


@pytest.fixture
def groups():
    return {
        Group(
            {Atom("C", (0, 0, 0), 0), Atom("O", (0, 0, 1), 1)},
            ConnectivityTable({(0, 1)}),
        ),
        Group(
            {Atom("O", (0, 1, 0), 2), Atom("H", (0, 1, 1), 3)},
            ConnectivityTable({(2, 3)}),
        ),
        Group(
            {Atom("O", (1, 0, 0), 4), Atom("H", (1, 0, 1), 5)},
            ConnectivityTable({(4, 5)}),
        ),
    }


@pytest.fixture
def bonds():
    return ConnectivityTable({(0, 1), (0, 2), (0, 4), (2, 3), (4, 5)})


def test_connect_groups(groups, bonds):
    assert connect_groups(groups, bonds) == ConnectivityTable({(0, 1), (0, 2)})
