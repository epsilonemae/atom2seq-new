import pytest

from atom2seq.atom_class import Atom
from atom2seq.cluster_class import Cluster
from atom2seq.fxnal_group_class import FxnalGroup


@pytest.fixture
def carbonyl():
    reset = FxnalGroup(set([]), set([]))
    reset._reset_used_indices()
    carbonyl = FxnalGroup(
        {
            Cluster({Atom("C", (0, 0, 0))}, set([])),
            Cluster({Atom("O", (1, 0, 0))}, set([])),
        },  # noqa
        {(1, 3)},
    )
    return carbonyl


@pytest.fixture
def carboxylic_acid():
    reset = FxnalGroup(set([]), set([]))
    reset._reset_used_indices()
    carboxylic_acid = FxnalGroup(
        {
            Cluster({Atom("C", (0, 0, 1))}, set([])),
            Cluster({Atom("O", (1, 0, 1))}, set([])),
            Cluster({Atom("O", (0, 1, 1)), Atom("H", (0, 1, 1))}, {(4, 5)}),
        },
        {(1, 3), (1, 6)},
    )
    return carboxylic_acid


@pytest.fixture
def amide():
    reset = FxnalGroup(set([]), set([]))
    reset._reset_used_indices()
    amide = FxnalGroup(
        {
            Cluster({Atom("C", (0, 0, 2))}, set([])),
            Cluster({Atom("O", (1, 0, 2))}, set([])),
            Cluster(
                {
                    Atom("N", (0, 1, 2)),
                    Atom("H", (1, 2, 2)),
                    Atom("H", (1, 1, 2)),
                },  # noqa
                {(4, 5), (4, 6)},
            ),
        },
        {(1, 3), (1, 7)},
    )
    return amide


@pytest.fixture
def phenyl():
    reset = FxnalGroup(set([]), set([]))
    reset._reset_used_indices()
    phenyl = FxnalGroup(
        {
            Cluster({Atom("C", (0, 0, 3)), Atom("H", (0, 1, 3))}, {(0, 1)}),
            Cluster({Atom("C", (1, 0, 3)), Atom("H", (1, 1, 3))}, {(3, 4)}),
            Cluster({Atom("C", (2, 0, 3)), Atom("H", (2, 1, 3))}, {(6, 7)}),
            Cluster({Atom("C", (3, 0, 3)), Atom("H", (3, 1, 3))}, {(9, 10)}),
            Cluster({Atom("C", (4, 0, 3)), Atom("H", (4, 1, 3))}, {(12, 13)}),
            Cluster({Atom("C", (5, 0, 3))}, set([])),
        },
        {(2, 5), (5, 8), (8, 11), (11, 14), (14, 16), (2, 16)},
    )
    return phenyl


@pytest.fixture
def phenol():
    reset = FxnalGroup(set([]), set([]))
    reset._reset_used_indices()
    phenol = FxnalGroup(
        {
            Cluster({Atom("C", (0, 0, 4)), Atom("H", (0, 1, 4))}, {(0, 1)}),
            Cluster({Atom("C", (1, 0, 4)), Atom("H", (1, 1, 4))}, {(3, 4)}),
            Cluster({Atom("C", (2, 0, 4))}, set([])),
            Cluster({Atom("C", (3, 0, 4)), Atom("H", (3, 1, 4))}, {(8, 9)}),
            Cluster({Atom("C", (4, 0, 4)), Atom("H", (4, 1, 4))}, {(11, 12)}),
            Cluster({Atom("C", (5, 0, 4))}, set([])),
            Cluster({Atom("O", (5, 1, 4)), Atom("H", (5, 2, 4))}, {(16, 17)}),
        },
        {(2, 5), (5, 7), (7, 10), (10, 13), (13, 15), (2, 15), (15, 18)},
    )
    return phenol


def test_get_symbol(carbonyl, carboxylic_acid, amide, phenyl, phenol):
    print(f"{carbonyl.get_symbol()=}")
    carbonyl._print_symbol()
    assert (
        (carbonyl.get_symbol() == "C=O")
        and (carboxylic_acid.get_symbol() == "COOH")
        and (amide.get_symbol() == "Amd")
        and (phenyl.get_symbol() == "Ph")
        and (phenol.get_symbol() == "PhOH")
    )


def test_rep_assignment(carbonyl, carboxylic_acid, amide, phenyl, phenol):
    assert (
        (carbonyl.get_rep() == 1)
        and (carboxylic_acid.get_rep() == 1)
        and (amide.get_rep() == 1)
        and (phenyl.get_rep() == 16)
        and (phenol.get_rep() == 7)
    )


def test_get_atoms(carbonyl):
    assert carbonyl.get_atoms() == {Atom("C", (0, 0, 0)), Atom("O", (1, 0, 0))}


def test_merge_clusters(carbonyl):
    carbonyl.merge_clusters(1, 3)
    # instead of dealing with the complications of equality with unique
    # indices, we check the repr to check it got merged properly.
    assert (
        repr(carbonyl)
        == "Group({Cluster({Atom('O', (1, 0, 0), 5, 2), Atom('C', (0, 0, 0), "
        "5, 0)}, {(0, 2)}, 4, 5)}, set(), -1, 4)"
    )


def test_dist(carbonyl, carboxylic_acid):
    assert carbonyl.dist(carboxylic_acid) == 1
