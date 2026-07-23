import pytest

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.get_pseq import Grouped, get_pseq
from atom2seq.group_class import Group


@pytest.fixture
def gap_groups():
    return {
        Group(
            {Atom("N", (0, 0, 0)), Atom("H", (0, 0, 1)), Atom("H", (0, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),  # idx 4
        Group(
            {Atom("C", (1, 0, 0)), Atom("H", (1, 0, 1)), Atom("H", (1, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),  # idx 5
        Group(
            {
                Atom("C", (2, 0, 0)),
                Atom("O", (2, 0, 1)),
                Atom("N", (2, 1, 0)),
                Atom("H", (2, 1, 1)),
            },
            ConnectivityTable({(0, 1), (0, 2), (2, 3)}),
        ),  # idx 9
        Group(
            {Atom("C", (3, 0, 0)), Atom("H", (3, 0, 1))},
            ConnectivityTable({(0, 1)}),
        ),  # idx 1
        Group(
            {
                Atom("C", (4, 0, 0)),
                Atom("H", (4, 0, 1)),
                Atom("H", (4, 1, 0)),
                Atom("H", (4, 1, 1)),
            },
            ConnectivityTable({(0, 1), (0, 2), (0, 3)}),
        ),  # idx 10
        Group(
            {Atom("C", (5, 0, 0)), Atom("O", (5, 0, 1))},
            ConnectivityTable({(0, 1)}),
        ),  # idx 2
        Group(
            {Atom("N", (6, 0, 0))},
            ConnectivityTable(set()),
        ),  # idx 0
        Group(
            {Atom("C", (7, 0, 0)), Atom("H", (7, 0, 1)), Atom("H", (7, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),  # idx 6
        Group(
            {Atom("C", (8, 0, 0)), Atom("H", (8, 0, 1)), Atom("H", (8, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),  # idx 7
        Group(
            {Atom("C", (9, 0, 0)), Atom("H", (9, 0, 1)), Atom("H", (9, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),  # idx 8
        Group(
            {Atom("C", (10, 0, 0)), Atom("H", (10, 0, 1))},
            ConnectivityTable({(0, 1)}),
        ),  # idx 3
        Group(
            {
                Atom("C", (11, 0, 0)),
                Atom("O", (11, 0, 1)),
                Atom("O", (11, 1, 0)),
                Atom("H", (11, 1, 1)),
            },
            ConnectivityTable({(0, 1), (0, 2), (2, 3)}),
        ),  # idx 11
    }


@pytest.fixture
def gap_bonds():
    return ConnectivityTable(
        {
            (0, 2),
            (0, 3),
            (0, 6),
            (1, 2),
            (1, 9),
            (1, 10),
            (3, 8),
            (3, 11),
            (4, 5),
            (5, 9),
            (6, 7),
            (7, 8),
        }
    )


@pytest.fixture
def gap(gap_groups, gap_bonds):
    return Grouped(gap_groups, gap_bonds)


@pytest.fixture
def no_nterm_groups():
    return {
        Group(
            {Atom("P", (0, 0, 0)), Atom("H", (0, 0, 1)), Atom("H", (0, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),
        Group(
            {Atom("C", (1, 0, 0)), Atom("H", (1, 0, 1)), Atom("H", (1, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),
        Group(
            {
                Atom("C", (2, 0, 0)),
                Atom("O", (2, 0, 1)),
                Atom("O", (2, 1, 0)),
                Atom("H", (2, 1, 1)),
            },
            ConnectivityTable({(0, 1), (0, 2), (2, 3)}),
        ),
    }


@pytest.fixture
def no_nterm_bonds():
    return ConnectivityTable({(0, 1), (1, 2)})


@pytest.fixture
def invalid_groups():
    return {
        Group(
            {Atom("N", (0, 0, 0)), Atom("H", (0, 0, 1)), Atom("H", (0, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),
        Group(
            {Atom("C", (1, 0, 0)), Atom("H", (1, 0, 1)), Atom("H", (1, 1, 0))},
            ConnectivityTable({(0, 1), (0, 2)}),
        ),
        Group(
            {
                Atom("C", (2, 0, 0)),
                Atom("S", (2, 0, 1)),
                Atom("O", (2, 1, 0)),
                Atom("H", (2, 1, 1)),
            },
            ConnectivityTable({(0, 1), (0, 2), (2, 3)}),
        ),
    }


@pytest.fixture
def invalid_bonds():
    return ConnectivityTable({(0, 1), (1, 2)})


def test_group_symbols(gap):
    assert gap.group_symbols(4) == ["H", "H", "N"]


def test_id_nterms(gap):
    assert gap.id_nterms() == [4]


def test_id_nterms_proline():
    nh = Grouped(
        {
            Group(
                {Atom("N", (0, 0, 0)), Atom("H", (0, 0, 1))},
                ConnectivityTable({(0, 1)}),
            )
        },
        ConnectivityTable(set()),
    )
    assert nh.id_nterms() == [0]


# Also tests rgroup_ider
def test_bb_iter(gap):
    assert gap.bb_iter(4) == ["G", "A", "P"]


def test_bb_iter_none(gap):
    assert gap.bb_iter(1) is None


def test_get_pseq(gap_groups, gap_bonds):
    assert get_pseq(gap_groups, gap_bonds) == ["G", "A", "P"]


def test_get_pseq_no_nterm_err(no_nterm_groups, no_nterm_bonds):
    with pytest.raises(ValueError) as err:
        get_pseq(no_nterm_groups, no_nterm_bonds)
    assert (
        err.exconly()
        == "ValueError: This is not a valid protein, as it has no potential "
        "N-termini."
    )


def test_get_pseq_no_backbone_err(invalid_groups, invalid_bonds):
    with pytest.raises(ValueError) as err:
        get_pseq(invalid_groups, invalid_bonds)
    assert (
        err.exconly()
        == "ValueError: This is not a valid protein, as the backbone could "
        "not be found."
    )
