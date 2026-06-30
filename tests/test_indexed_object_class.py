import pytest

from atom2seq.indexed_object_class import IndexedObject as IO


@pytest.fixture
def ios():
    IO()._reset_used_indices()
    io0 = IO()
    io1 = IO()
    io2 = IO()
    io3 = IO()
    return (io0, io1, io2, io3)


def test_eq(ios):
    assert (ios[0] == ios[0]) and (ios[0] != ios[1])


def test_lt(ios):
    assert (ios[0] < ios[1]) and (ios[2] > ios[1])


def test_tuple(ios):
    assert ios[0]._tuple() == (0,)


def test_idx_neq(ios):
    assert not ios[0]._idx_neq(ios[1])


def test_get_idx(ios):
    assert (ios[0].get_idx() == 0) and (ios[1].get_idx() == 1)


def test_set_idx(ios):
    ios[0].set_idx(4)
    assert ios[0].get_idx() == 4
