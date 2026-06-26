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


def test_get_atoms(clusters):
    cluster1, cluster2 = clusters
    assert cluster1.get_atoms() == {
        Atom("C", (0, 0, 0)),
        Atom("H", (1, 0, 0)),
        Atom("H", (0, 1, 0)),
    }


def test_eq(clusters):
    cluster1, cluster2 = clusters
    assert (
        (cluster1 == cluster1)
        and (cluster2 == cluster2)
        and (cluster1 != cluster2)  # noqa
    )  # noqa


def test_set_atoms(clusters):
    cluster1, cluster2 = clusters
    cluster1.set_atoms({Atom("S", (0, 0, 2)), Atom("H", (1, 0, 2))})
    assert (cluster1.get_atoms() == cluster2.get_atoms()) and (
        not cluster1.get_bonds()
    )  # noqa


def test_add_atom(clusters):
    cluster1, cluster2 = clusters
    cluster1.add_atom("H", (0, 0, 1))
    assert (
        cluster1.get_atoms()
        == {
            Atom("C", (0, 0, 0)),
            Atom("H", (1, 0, 0)),
            Atom("H", (0, 1, 0)),
            Atom("H", (0, 0, 1)),
        }
    ) and (cluster1.get_symbol() == "CH3")


def test_del_atom(clusters):
    cluster1, cluster2 = clusters
    cluster1.del_atom(2)
    assert (
        cluster1.get_atoms()
        == {
            Atom("C", (0, 0, 0)),
            Atom("H", (1, 0, 0)),
        }
    ) and cluster1.get_symbol() == "CH"


def test_get_bonds(clusters):
    cluster1, cluster2 = clusters
    assert cluster1.get_bonds() == {(0, 1), (0, 2)}


def test_set_bonds(clusters):
    cluster1, cluster2 = clusters
    cluster1.set_bonds({(0, 2), (1, 2)})
    assert cluster1.get_bonds() == {(0, 2), (1, 2)}


def test_add_bond(clusters):
    cluster1, cluster2 = clusters
    cluster1.add_bond(1, 2)
    assert cluster1.get_bonds() == {(0, 1), (0, 2), (1, 2)}


def test_del_bond(clusters):
    cluster1, cluster2 = clusters
    cluster1.del_bond(0, 2)
    assert cluster1.get_bonds() == {(0, 1)}


def test_check_bond(clusters):
    cluster1, cluster2 = clusters
    assert (cluster1.check_bond(0, 1)) and (not cluster1.check_bond(1, 2))


def test_get_parent(clusters):
    cluster1, cluster2 = clusters
    assert cluster1.get_parent() == cluster2.get_parent() == -1


def test_set_parent(clusters):
    cluster1, cluster2 = clusters
    cluster1.set_parent(1)
    cluster2.set_parent(2)
    assert (cluster1.get_parent() == 1) and (cluster2.get_parent() == 2)


def test_get_rep(clusters):
    cluster1, cluster2 = clusters
    assert (cluster1.get_rep() == cluster1.used_indices[0]) and (
        cluster2.get_rep() == cluster2.used_indices[4]
    )


def test_set_rep(clusters):
    cluster1, cluster2 = clusters
    print(cluster1.used_indices)
    print(list(cluster1.get_atoms())[0].used_indices)
    cluster1.set_rep(1)
    cluster2.set_rep(5)
    assert (cluster1.get_rep() == cluster1.used_indices[1]) and (
        cluster2.get_rep() == cluster2.used_indices[5]
    )


def test_dist(clusters):
    cluster1, cluster2 = clusters
    assert cluster1.dist(cluster2) == 2
