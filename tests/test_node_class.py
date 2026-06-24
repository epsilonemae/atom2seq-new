import pytest

from atom2seq.classes import Atom
from atom2seq.classes import IndexedObject as IO
from atom2seq.classes import Node


@pytest.fixture
def nodes():
    IO.used_indices = {}
    node1 = Node(
        {
            Atom("C", (0, 0, 0)),
            Atom("H", (1, 0, 0)),
            Atom("H", (0, 1, 0)),
        },
        {(0, 1), (0, 2)},
    )
    node2 = Node({Atom("S", (0, 0, 2)), Atom("H", (1, 0, 2))}, {(4, 5)})
    return (node1, node2)


def test_get_symbol(nodes):
    node1, node2 = nodes
    assert (node1.get_symbol() == "CH2") and (node2.get_symbol() == "SH")


def test_get_atoms(nodes):
    node1, node2 = nodes
    assert node1.get_atoms() == {
        Atom("C", (0, 0, 0)),
        Atom("H", (1, 0, 0)),
        Atom("H", (0, 1, 0)),
    }


def test_eq(nodes):
    node1, node2 = nodes
    assert (node1 == node1) and (node2 == node2) and (node1 != node2)


def test_set_atoms(nodes):
    node1, node2 = nodes
    node1.set_atoms({Atom("S", (0, 0, 2)), Atom("H", (1, 0, 2))})
    assert (node1.get_atoms() == node2.get_atoms()) and (not node1.get_bonds())


def test_add_atom(nodes):
    node1, node2 = nodes
    print(f"{IO.used_indices=}, {Atom.used_indices=}, {Node.used_indices=}")
    node1.add_atom("H", (0, 0, 1))
    assert (
        node1.get_atoms()
        == {
            Atom("C", (0, 0, 0)),
            Atom("H", (1, 0, 0)),
            Atom("H", (0, 1, 0)),
            Atom("H", (0, 0, 1)),
        }
    ) and (node1.get_symbol() == "CH3")


def test_del_atom(nodes):
    node1, node2 = nodes
    node1.del_atom(2)
    assert (
        node1.get_atoms()
        == {
            Atom("C", (0, 0, 0)),
            Atom("H", (1, 0, 0)),
        }
    ) and node1.get_symbol() == "CH"


def test_get_bonds(nodes):
    node1, node2 = nodes
    assert node1.get_bonds() == {(0, 1), (0, 2)}


def test_set_bonds(nodes):
    node1, node2 = nodes
    node1.set_bonds({(0, 2), (1, 2)})
    assert node1.get_bonds() == {(0, 2), (1, 2)}


def test_add_bond(nodes):
    node1, node2 = nodes
    node1.add_bond(1, 2)
    assert node1.get_bonds() == {(0, 1), (0, 2), (1, 2)}


def test_del_bond(nodes):
    node1, node2 = nodes
    node1.del_bond(0, 2)
    assert node1.get_bonds() == {(0, 1)}


def test_check_bond(nodes):
    node1, node2 = nodes
    assert (node1.check_bond(0, 1)) and (not node1.check_bond(1, 2))


def test_get_idx(nodes):
    node1, node2 = nodes
    assert (node1.get_idx() == 3) and (node2.get_idx() == 6)


def test_set_idx(nodes):
    node1, node2 = nodes
    node1.set_idx(7)
    node2.set_idx(8)
    assert (node1.get_idx() == 7) and (node2.get_idx() == 8)


def test_get_parent(nodes):
    node1, node2 = nodes
    assert node1.get_parent() == node2.get_parent() == -1


def test_set_parent(nodes):
    node1, node2 = nodes
    node1.set_parent(1)
    node2.set_parent(2)
    assert (node1.get_parent() == 1) and (node2.get_parent() == 2)


def test_get_rep(nodes):
    node1, node2 = nodes
    assert (node1.get_rep() == node1.used_indices[0]) and (
        node2.get_rep() == node2.used_indices[4]
    )


def test_set_rep(nodes):
    node1, node2 = nodes
    node1.set_rep(1)
    node2.set_rep(5)
    assert (node1.get_rep() == node1.used_indices[1]) and (
        node2.get_rep() == node2.used_indices[5]
    )


def test_dist(nodes):
    node1, node2 = nodes
    assert node1.dist(node2) == 2
