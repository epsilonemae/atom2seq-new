import pytest

from atom2seq.backbone_class import Backbone
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


@pytest.fixture
def empty_bb():
    return Backbone(
        {Group(set(), ConnectivityTable(set()))}, ConnectivityTable(set())
    )  # noqa


def test_repr(empty_bb):
    assert (
        repr(empty_bb) == "Backbone({Group([], ConnectivityTable([]))}, "
        "ConnectivityTable([]))"
    )


def test_eq(empty_bb):
    assert empty_bb == Backbone(
        {Group(set(), ConnectivityTable(set()))}, ConnectivityTable(set())
    )


def test_get_bonds(empty_bb):
    assert empty_bb.get_bonds() == ConnectivityTable(set())


def test_get_groups(empty_bb):
    assert empty_bb.get_groups() == {Group(set(), ConnectivityTable(set()))}


def test_group_list(empty_bb):
    assert empty_bb.group_list() == [Group(set(), ConnectivityTable(set()))]


def test_add_group(empty_bb):
    empty_bb.add_group(Group(set(), ConnectivityTable(set())))
    assert empty_bb.get_groups() == {
        Group(set(), ConnectivityTable(set())),
        Group(set(), ConnectivityTable(set())),
    }


def test_add_group_warning(empty_bb):
    with pytest.warns(UserWarning):
        group_to_add = Group(set(), ConnectivityTable(set()))
        group_to_add.set_idx(0)
        empty_bb.add_group(group_to_add)
        assert empty_bb.get_groups() == {
            Group(set(), ConnectivityTable(set())),
            Group(set(), ConnectivityTable(set())),
        }
