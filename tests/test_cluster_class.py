import pytest

from atom2seq.atom_class import Atom
from atom2seq.cluster_class import Cluster


@pytest.fixture
def clusters():
    reset = Cluster(set([]), set([]))
    reset._reset_used_indices()
    cluster1 = Cluster(
        {
            Atom("C", (0, 0, 0)),
            Atom("H", (1, 0, 0)),
            Atom("H", (0, 1, 0)),
        },
        {(0, 1), (0, 2)},
    )
    cluster2 = Cluster({Atom("S", (0, 0, 2)), Atom("H", (1, 0, 2))}, {(4, 5)})
    return (cluster1, cluster2)


def test_get_symbol(clusters):
    cluster1, cluster2 = clusters
    assert (cluster1.get_symbol() == "CH2") and (cluster2.get_symbol() == "SH")


def test_dist(clusters):
    cluster1, cluster2 = clusters
    assert cluster1.dist(cluster2) == 2
