import pytest

from atom2seq.connectivity_table_class import ConnectivityTable


@pytest.fixture
def graph():
    return ConnectivityTable({(0, 1), (1, 2), (2, 3), (2, 4)})


def test_eq(graph):
    assert graph == ConnectivityTable({(0, 1), (1, 2), (2, 3), (2, 4)})


def test_len(graph):
    assert len(graph) == 4


def test_get_pairs(graph):
    assert graph.get_pairs() == {(0, 1), (1, 2), (2, 3), (2, 4)}


def test_add_pair(graph):
    graph.add_pair((3, 5))
    assert graph.get_pairs() == {(0, 1), (1, 2), (2, 3), (2, 4), (3, 5)}


def test_del_pair(graph):
    graph.del_pair((0, 1))
    assert graph.get_pairs() == {(1, 2), (2, 3), (2, 4)}


def test_check_pair(graph):
    assert (graph.check_pair((0, 1))) and (not graph.check_pair((3, 4)))


def test_get_paired(graph):
    assert graph.get_paired(2) == {1, 3, 4}
