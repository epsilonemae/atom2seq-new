import pytest

from atom2seq.atom_class import Atom
from atom2seq.cluster_class import Cluster
from atom2seq.group_class import Group


@pytest.fixture
def carbonyl():
    reset = Group(set([]), set([]))
    reset._reset_used_indices()
    carbonyl = Group(
        {
            Cluster({Atom("C", (0, 0, 0))}, set([])),
            Cluster({Atom("O", (1, 0, 0))}, set([])),
        },  # noqa
        {(1, 3)},
    )
    return carbonyl


@pytest.fixture
def carboxylic_acid():
    reset = Group(set([]), set([]))
    reset._reset_used_indices()
    carboxylic_acid = Group(
        {
            Cluster({Atom("C", (0, 0, 0))}, set([])),
            Cluster({Atom("O", (1, 0, 0))}, set([])),
            Cluster({Atom("O", (0, 1, 0)), Atom("H", (0, 1, 1))}, {(4, 5)}),
        },
        {(1, 3), (1, 6)},
    )
    return carboxylic_acid


@pytest.fixture
def amide():
    reset = Group(set([]), set([]))
    reset._reset_used_indices()
    amide = Group(
        {
            Cluster({Atom("C", (0, 0, 0))}, set([])),
            Cluster({Atom("O", (1, 0, 0))}, set([])),
            Cluster(
                {
                    Atom("N", (0, 1, 0)),
                    Atom("H", (0, 1, 1)),
                    Atom("H", (1, 1, 0)),
                },  # noqa
                {(4, 5), (4, 6)},
            ),
        },
        {(1, 3), (1, 7)},
    )
    return amide


@pytest.fixture
def phenyl():
    reset = Group(set([]), set([]))
    reset._reset_used_indices()
    phenyl = Group(
        {
            Cluster({Atom("C", (0, 0, 0)), Atom("H", (0, 1, 0))}, {(0, 1)}),
            Cluster({Atom("C", (1, 0, 0)), Atom("H", (1, 1, 0))}, {(3, 4)}),
            Cluster({Atom("C", (2, 0, 0)), Atom("H", (2, 1, 0))}, {(6, 7)}),
            Cluster({Atom("C", (3, 0, 0)), Atom("H", (3, 1, 0))}, {(9, 10)}),
            Cluster({Atom("C", (4, 0, 0)), Atom("H", (4, 1, 0))}, {(12, 13)}),
            Cluster({Atom("C", (5, 0, 0))}, set([])),
        },
        {(2, 5), (5, 8), (8, 11), (11, 14), (14, 16), (2, 16)},
    )
    return phenyl


@pytest.fixture
def phenol():
    reset = Group(set([]), set([]))
    reset._reset_used_indices()
    phenol = Group(
        {
            Cluster({Atom("C", (0, 0, 0)), Atom("H", (0, 1, 0))}, {(0, 1)}),
            Cluster({Atom("C", (1, 0, 0)), Atom("H", (1, 1, 0))}, {(3, 4)}),
            Cluster({Atom("C", (2, 0, 0))}, set([])),
            Cluster({Atom("C", (3, 0, 0)), Atom("H", (3, 1, 0))}, {(8, 9)}),
            Cluster({Atom("C", (4, 0, 0)), Atom("H", (4, 1, 0))}, {(11, 12)}),
            Cluster({Atom("C", (5, 0, 0))}, set([])),
            Cluster({Atom("O", (5, 0, 1)), Atom("H", (5, 1, 1))}, {(16, 17)}),
        },
        {(2, 5), (5, 7), (7, 10), (10, 13), (13, 15), (2, 15), (15, 18)},
    )
    return phenol


def test_get_symbol(carbonyl, carboxylic_acid, amide, phenyl, phenol):
    assert (
        (carbonyl.get_symbol() == "C=O")
        and (carboxylic_acid.get_symbol() == "COOH")
        and (amide.get_symbol() == "Amd")
        and (phenyl.get_symbol() == "Ph")
        and (phenol.get_symbol() == "PhOH")
    )


def test_get_clusters(carbonyl):
    assert carbonyl.get_clusters() == {
        Cluster({Atom("C", (0, 0, 0))}, set([])),
        Cluster({Atom("O", (1, 0, 0))}, set([])),
    }


def test_set_clusters(carbonyl):
    carbonyl.set_clusters(
        {
            Cluster({Atom("C", (1, 0, 0))}, set([])),
            Cluster({Atom("O", (2, 0, 0))}, set([])),
        }
    )
    assert (
        carbonyl.get_clusters()
        == {
            Cluster({Atom("C", (1, 0, 0))}, set([])),
            Cluster({Atom("O", (2, 0, 0))}, set([])),
        }
    ) and (not carbonyl.get_bonds())


def test_add_cluster(carbonyl):
    carbonyl.add_cluster({Atom(("N", (1, 1, 1)))}, set([]))
    assert carbonyl.get_clusters() == {
        Cluster({Atom("C", (0, 0, 0))}, set([])),
        Cluster({Atom("O", (1, 0, 0))}, set([])),
        Cluster({Atom("N", (1, 1, 1))}, set([])),
    }


def test_del_cluster(carbonyl):
    carbonyl.del_cluster(1)
    assert carbonyl.get_clusters() == {
        Cluster({Atom("O", (1, 0, 0))}, set([])),
    }


def test_merge_clusters(carbonyl):
    carbonyl.merge_clusters(1, 3)
    assert carbonyl.get_clusters() == {
        Cluster({Atom("C", (1, 0, 0)), Atom("O", (2, 0, 0))}, {0, 1})
    }


def test_get_bonds(carbonyl):
    assert carbonyl.get_bonds() == {(1, 3)}


def test_set_bonds(carboxylic_acid):
    carboxylic_acid.set_bonds({(1, 3), (3, 6)})
    assert carboxylic_acid.get_bonds() == {(1, 3), (3, 6)}


def test_add_bond(carboxylic_acid):
    carboxylic_acid.add_bond(3, 6)
    assert carboxylic_acid.get_bonds() == {(1, 3), (1, 6), (3, 6)}


def test_del_bond(carbonyl):
    carbonyl.del_bond(1, 3)
    assert carbonyl.get_bonds() == set([])


def test_check_bond(carbonyl):
    assert carbonyl.check_bond(1, 3)


def get_rep(carbonyl):
    assert carbonyl.get_rep() == 1


def set_rep(carboxylic_acid):
    carboxylic_acid.set_rep(6)
    assert carboxylic_acid.get_rep() == 6
