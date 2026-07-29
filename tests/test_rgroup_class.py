import pytest

from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.rgroup_class import RGroup


@pytest.fixture
def empty_rgroup():
    return RGroup(
        {Group(set(), ConnectivityTable(set()))}, ConnectivityTable(set())
    )  # noqa


def test_repr(empty_rgroup):
    assert (
        repr(empty_rgroup) == "RGroup([Group([], ConnectivityTable([]), 0)], "
        "ConnectivityTable([]))"
    )


def test_eq(empty_rgroup):
    assert empty_rgroup == RGroup(
        {Group(set(), ConnectivityTable(set()))}, ConnectivityTable(set())
    )


def test_get_bonds(empty_rgroup):
    assert empty_rgroup.get_bonds() == ConnectivityTable(set())


def test_get_groups(empty_rgroup):
    assert empty_rgroup.get_groups() == {Group(set(), ConnectivityTable(set()))}  # noqa


def test_group_list(empty_rgroup):
    assert empty_rgroup.group_list() == [Group(set(), ConnectivityTable(set()))]  # noqa


def test_add_group(empty_rgroup):
    empty_rgroup.add_group(Group(set(), ConnectivityTable(set())))
    assert empty_rgroup.get_groups() == {
        Group(set(), ConnectivityTable(set())),
        Group(set(), ConnectivityTable(set())),
    }
