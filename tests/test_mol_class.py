# import pytest

# from atom2seq.atom_class import Atom
# from atom2seq.cluster_class import Cluster
# from atom2seq.fxnal_group_class import FxnalGroup
# from atom2seq.mol_class import Mol


# @pytest.fixture
# def dih2():
#     molecule = Mol(
#         {
#             FxnalGroup({Cluster({Atom("H", (0, 0, 0))}, {})}, {}),
#             FxnalGroup({Cluster({Atom("H", (0, 0, 1))}, {})}, {}),
#             FxnalGroup({Cluster({Atom("H", (0, 0, 2))}, {})}, {}),
#             FxnalGroup({Cluster({Atom("H", (0, 0, 3))}, {})}, {}),
#         },
#         {(2, 5), (8, 11)},
#     )
#     return molecule


# def test_get_backbone(dih2):
#     assert dih2.get_backbone() == set()


# def test_set_backbone(dih2):
#     dih2.set_backbone({2, 5})
#     assert dih2.get_backbone() == {2, 5}


# def test_add_backbone(dih2):
#     dih2.add_backbone(2)
#     assert dih2.get_backbone() == {
#         2,
#     }


# def test_del_backbone(dih2):
#     dih2.set_backbone({2, 5})
#     dih2.del_backbone(5)
#     assert dih2.get_backbone() == {
#         2,
#     }


# def test_check_backbone(dih2):
#     dih2.set_backbone({2, 5})
#     assert (dih2.check_backbone(2)) and (not dih2.check_backbone(8))


# def test_get_nterm(dih2):
#     assert dih2.get_nterm() == -1


# def test_set_nterm(dih2):
#     dih2.set_nterm(2)
#     assert dih2.get_nterm() == 2


# def test_get_atoms(dih2):
#     assert dih2.get_atoms() == {
#         Atom("H", (0, 0, 0)),
#         Atom("H", (0, 0, 1)),
#         Atom("H", (0, 0, 2)),
#         Atom("H", (0, 0, 3)),
#     }


# def test_get_clusters(dih2):
#     assert dih2.get_clusters() == {
#         Cluster({Atom("H", (0, 0, 0))}, {}),
#         Cluster({Atom("H", (0, 0, 1))}, {}),
#         Cluster({Atom("H", (0, 0, 2))}, {}),
#         Cluster({Atom("H", (0, 0, 3))}, {}),
#     }


# def test_find_submol(dih2):
#     assert dih2.find_submol(2) == {2, 5}


# def test_del_submol(dih2):
#     dih2.del_submol(2)
#     assert dih2.get_vertices() == {
#         FxnalGroup({Cluster({Atom("H", (0, 0, 2))}, {})}, {}),
#         FxnalGroup({Cluster({Atom("H", (0, 0, 3))}, {})}, {}),
#     }


# def test_pop_submol(dih2):
#     submol = dih2.pop_submol(2)
#     assert (
#         dih2.get_vertices()
#         == {
#             FxnalGroup({Cluster({Atom("H", (0, 0, 2))}, {})}, {}),
#             FxnalGroup({Cluster({Atom("H", (0, 0, 3))}, {})}, {}),
#         }
#     ) and (
#         submol.get_vertices()
#         == {
#             FxnalGroup({Cluster({Atom("H", (0, 0, 0))}, {})}, {}),
#             FxnalGroup({Cluster({Atom("H", (0, 0, 1))}, {})}, {}),
#         }
#     )
