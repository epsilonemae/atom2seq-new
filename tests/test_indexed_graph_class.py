import pytest

from atom2seq.indexed_graph_class import IndexedGraph
from atom2seq.indexed_object_class import IndexedObject as IO


@pytest.fixture
def graph():
    IO()._reset_used_indices()
    graph = IndexedGraph(
        {IO(), IO(), IO()},
        {(0, 1), (0, 2)},
    )
    return graph


def test_vertex_indices(graph):
    assert graph.vertex_indices() == {0, 1, 2}


def test_vertex_list(graph):
    IO()._reset_used_indices()
    assert graph.vertex_list() == [
        IO(idx=0),
        IO(idx=1),
        IO(idx=2),
    ]


def test_get_vertices(graph):
    IO()._reset_used_indices()
    assert graph.get_vertices() == {
        IO(idx=0),
        IO(idx=1),
        IO(idx=2),
    }


def test_set_vertices(graph):
    graph.set_vertices({IO()})
    IO()._reset_used_indices()
    assert (graph.get_vertices() == {IO(idx=3)}) and (not graph.get_edges())


def test_add_vertex(graph):
    graph.add_vertex(IO())
    IO()._reset_used_indices()
    assert graph.get_vertices() == {
        IO(idx=0),
        IO(idx=1),
        IO(idx=2),
        IO(idx=3),
    }


def test_del_vertex(graph):
    graph.del_vertex(2)
    IO()._reset_used_indices()
    assert graph.get_vertices() == {
        IO(idx=0),
        IO(idx=1),
    }


def test_get_edges(graph):
    assert graph.get_edges() == {(0, 1), (0, 2)}


def test_set_edges_valid(graph):
    graph.set_edges({(0, 1), (1, 2)})
    assert graph.get_edges() == {(0, 1), (1, 2)}


def test_set_edges_invalid(graph):
    with pytest.raises(ValueError) as error_msg:
        graph.set_edges({(0, 1), (1, 3)})
        assert (
            error_msg == "The edge (1, 3) contains the index of an object "
            "that is not a child of this graph."
        )


def test_add_edge(graph):
    graph.add_edge((1, 2))
    assert graph.get_edges() == {(0, 1), (0, 2), (1, 2)}


def test_del_edge(graph):
    graph.del_edge((0, 2))
    assert graph.get_edges() == {(0, 1)}


def test_check_edge(graph):
    assert graph.check_edge((0, 1)) and not graph.check_edge((1, 2))


def test_get_rep(graph):
    assert graph.get_rep() == -1


def test_set_rep_valid(graph):
    graph.set_rep(1)
    assert graph.get_rep() == 1


def test_set_rep_invalid(graph):
    with pytest.raises(ValueError) as error_msg:
        graph.set_rep(3)
        assert (
            error_msg == "The object at index 3 is not a child of this graph."
        )  # noqa


def test_get_adjacent(graph):
    assert graph.get_adjacent(0) == {1, 2}
