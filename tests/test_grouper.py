import pytest

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.grouper import GroupedMol, group_mol
from atom2seq.mol_class import Mol


@pytest.fixture
def mol():
    return Mol(
        {
            Atom("H", (0, 1, 0)),
            Atom("H", (1, 2, 0)),
            Atom("C", (1, 1, 0)),
            Atom("H", (1, 0, 0)),
            Atom("O", (2, 2, 0)),
            Atom("C", (2, 1, 0)),
            Atom("N", (3, 1, 0)),
            Atom("H", (3, 0, 0)),
            Atom("H", (4, 2, 0)),
            Atom("C", (4, 1, 0)),
            Atom("H", (4, 0, 0)),
            Atom("C", (5, 1, 0)),
            Atom("O", (5, 0, 0)),
            Atom("O", (6, 2, 0)),
            Atom("C", (6, 1, 0)),
            Atom("O", (7, 1, 0)),
            Atom("H", (8, 1, 0)),
        },
        ConnectivityTable(
            {
                (0, 2),
                (1, 2),
                (2, 3),
                (2, 4),
                (4, 5),
                (4, 7),
                (6, 7),
                (7, 9),
                (8, 9),
                (9, 10),
                (9, 12),
                (11, 12),
                (12, 13),
                (13, 14),
                (13, 15),
                (15, 16),
            }
        ),
    )


@pytest.fixture
def gmol(mol):
    return GroupedMol(mol)


def test_detectXHn(gmol):
    gmol.detectXHn(0)
    assert gmol.groups == {
        Group(
            {
                Atom("H", (0, 1, 0)),
                Atom("H", (1, 0, 0)),
                Atom("C", (1, 1, 0)),
                Atom("H", (1, 2, 0)),
            },
            ConnectivityTable({(0, 2), (1, 2), (2, 3)}),
        )
    }


def test_detectXHn_onlyH():
    test = GroupedMol(
        Mol(
            {
                Atom("C", (0, 0, 0)),
                Atom("O", (0, 0, 1)),
                Atom("H", (0, 1, 0)),
                Atom("H", (1, 0, 0)),
            },
            ConnectivityTable({(0, 1), (0, 2), (0, 3)}),
        )
    )
    test.detectCO(0)
    test.detectXHn(2)
    assert test.groups == {
        Group(
            {
                Atom("C", (0, 0, 0)),
                Atom("O", (0, 0, 1)),
            },
            ConnectivityTable({(0, 1)}),
        ),
        Group(
            {
                Atom("H", (0, 1, 0)),
            },
            ConnectivityTable(set()),
        ),
    }


def test_detectCO(gmol):
    gmol.detectCO(5)
    gmol.detectCO(12)
    gmol.detectCO(14)
    assert gmol.groups == {
        Group(
            {Atom("C", (2, 1, 0)), Atom("O", (2, 2, 0))},
            ConnectivityTable({(4, 5)}),  # noqa
        ),
        Group(
            {Atom("C", (5, 1, 0)), Atom("O", (5, 0, 0))},
            ConnectivityTable({(11, 12)}),  # noqa
        ),
        Group(
            {Atom("C", (6, 1, 0)), Atom("O", (6, 2, 0))},
            ConnectivityTable({(13, 14)}),  # noqa
        ),
    }


def test_detect_CO_notCO(gmol):
    gmol.detectCO(9)  # this is a C in a CH2
    assert gmol.groups == set()


def test_detectCOOH(gmol):
    gmol.detectCOOH(16)  # run it on the H since it runs everything (H-O-C=O)
    assert gmol.groups == {
        Group(
            {
                Atom("O", (6, 2, 0)),
                Atom("C", (6, 1, 0)),
                Atom("O", (7, 1, 0)),
                Atom("H", (8, 1, 0)),
            },
            ConnectivityTable({(13, 14), (13, 15), (15, 16)}),
        )
    }


def test_detectCOOH_notCOOH(gmol):
    gmol.detectCOOH(4)  # this C is in an Amd not a COOH
    gmol.detectCOOH(5)  # this O is in an Amd not a COOH
    gmol.detectCOOH(6)  # this H is in an Amd not a COOH
    assert gmol.groups == set()


def test_detectAmd(gmol):
    gmol.detectAmd(6)  # running it on the H since it runs everything (H-N-C=O)
    assert gmol.groups == {
        Group(
            {
                Atom("O", (2, 2, 0)),
                Atom("C", (2, 1, 0)),
                Atom("N", (3, 1, 0)),
                Atom("H", (3, 0, 0)),
            },
            ConnectivityTable({(4, 5), (4, 7), (6, 7)}),
        )
    }


def test_detectAmd_notAmd(gmol):
    gmol.detectAmd(0)  # this H is in a CH3 not an Amd
    gmol.detectAmd(2)  # this C is in a CH3 not an Amd
    gmol.detectAmd(14)  # this O is in a COOH not an Amd
    assert gmol.groups == set()


def test_group_mol(mol):
    assert group_mol(mol) == {
        Group(
            {
                Atom("H", (0, 1, 0)),
                Atom("H", (1, 2, 0)),
                Atom("C", (1, 1, 0)),
                Atom("H", (1, 0, 0)),
            },
            ConnectivityTable({(0, 2), (1, 2), (2, 3)}),
        ),
        Group(
            {
                Atom("O", (2, 2, 0)),
                Atom("C", (2, 1, 0)),
                Atom("N", (3, 1, 0)),
                Atom("H", (3, 0, 0)),
            },
            ConnectivityTable(
                {
                    (4, 5),
                    (4, 7),
                    (6, 7),
                }
            ),
        ),
        Group(
            {
                Atom("H", (4, 2, 0)),
                Atom("C", (4, 1, 0)),
                Atom("H", (4, 0, 0)),
            },
            ConnectivityTable({(8, 9), (9, 10)}),
        ),
        Group(
            {
                Atom("C", (5, 1, 0)),
                Atom("O", (5, 0, 0)),
            },
            ConnectivityTable({(11, 12)}),
        ),
        Group(
            {
                Atom("O", (6, 2, 0)),
                Atom("C", (6, 1, 0)),
                Atom("O", (7, 1, 0)),
                Atom("H", (8, 1, 0)),
            },
            ConnectivityTable({(13, 14), (13, 15), (15, 16)}),
        ),
    }
