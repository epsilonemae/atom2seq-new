import pytest

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.rgroup_class import RGroup


# Coordinates and in-group connectivity are unimportant here so all atoms are
# going at (n, 0, 0) and all groups will have an empty ConnectivityTable. This
# will allow me to simply pass the dictionary key into this function and get
# out the set of groups, then make the ConnectivityTable for each AA.
def group_set(syms_lol):
    print(syms_lol)
    out = set()
    n = 0
    for sym_list in syms_lol:
        print(sym_list)
        atoms = set()
        for sym in sym_list:
            print(sym)
            atoms.add(Atom(sym, (n, 0, 0)))
            print(atoms)
            n += 1
        out.add(Group(atoms, ConnectivityTable(set())))
        print(out)
    return out


def test_A():
    ala = RGroup(
        group_set((("C", "H", "H", "H"),)),
        ConnectivityTable(set()),
    )
    assert ala.symbol() == "A"


def test_C():
    cys = RGroup(
        group_set((("C", "H", "H"), ("H", "S"))), ConnectivityTable({(0, 1)})
    )  # noqa
    assert cys.symbol() == "C"


def test_D():
    asp = RGroup(
        group_set((("C", "H", "H"), ("C", "H", "O", "O"))),
        ConnectivityTable({(0, 1)}),  # noqa
    )
    assert asp.symbol() == "D"


def test_E():
    glu = RGroup(
        group_set((("C", "H", "H"), ("C", "H", "H"), ("C", "H", "O", "O"))),
        ConnectivityTable({(0, 1), (1, 2)}),
    )
    assert glu.symbol() == "E"


def test_F():
    phe = RGroup(
        group_set(
            (
                ("C",),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H", "H"),
            )
        ),
        ConnectivityTable({(0, 1)}),
    )
    assert phe.symbol() == "F"


def test_H():
    his = RGroup(
        group_set(
            (
                ("C",),
                ("C", "H"),
                ("C", "H"),
                ("C", "H", "H"),
                ("H", "N"),
                ("N",),
            )
        ),
        ConnectivityTable(set()),
    )
    assert his.symbol() == "H"


def test_I():
    ile = RGroup(
        group_set(
            (
                ("C", "H"),
                ("C", "H", "H"),
                ("C", "H", "H", "H"),
                ("C", "H", "H", "H"),
            )
        ),
        ConnectivityTable({(0, 1), (0, 2), (2, 3)}),
    )
    assert ile.symbol() == "I"


def test_K():
    lys = RGroup(
        group_set(
            (
                ("C", "H", "H"),
                ("C", "H", "H"),
                ("C", "H", "H"),
                ("C", "H", "H"),
                ("H", "H", "N"),
            )
        ),
        ConnectivityTable({(0, 1), (1, 2), (2, 3), (3, 4)}),
    )
    assert lys.symbol() == "K"


def test_L():
    leu = RGroup(
        group_set(
            (
                ("C", "H"),
                ("C", "H", "H"),
                ("C", "H", "H", "H"),
                ("C", "H", "H", "H"),
            )
        ),
        ConnectivityTable({(0, 1), (0, 2), (0, 3)}),
    )
    assert leu.symbol() == "L"


def test_M():
    met = RGroup(
        group_set(
            (("C", "H", "H"), ("C", "H", "H"), ("C", "H", "H", "H"), ("S",))
        ),  # noqa
        ConnectivityTable({(0, 1), (1, 3), (2, 3)}),
    )
    assert met.symbol() == "M"


def test_N():
    asn = RGroup(
        group_set((("C", "H", "H"), ("C", "H", "N", "O"), ("H",))),
        ConnectivityTable({(0, 1)}),
    )
    assert asn.symbol() == "N"


def test_P():
    pro = RGroup(
        group_set((("C", "H", "H"), ("C", "H", "H"), ("C", "H", "H"))),
        ConnectivityTable({(0, 1), (1, 2)}),
    )
    assert pro.symbol() == "P"


def test_Q():
    gln = RGroup(
        group_set(
            (("C", "H", "H"), ("C", "H", "H"), ("C", "H", "N", "O"), ("H",))
        ),  # noqa
        ConnectivityTable({(0, 1), (1, 2)}),
    )
    assert gln.symbol() == "Q"


def test_R():
    arg = RGroup(
        group_set(
            (
                ("C",),
                ("C", "H", "H"),
                ("C", "H", "H"),
                ("C", "H", "H"),
                ("H", "H", "N"),
                ("H", "N"),
                ("H", "N"),
            )
        ),
        ConnectivityTable({(1, 2), (2, 3), (3, 5), (0, 5), (0, 4), (0, 6)}),
    )
    assert arg.symbol() == "R"


def test_S():
    ser = RGroup(
        group_set((("C", "H", "H"), ("H", "O"))), ConnectivityTable({(0, 1)})
    )  # noqa
    assert ser.symbol() == "S"


def test_T():
    thr = RGroup(
        group_set((("C", "H"), ("C", "H", "H", "H"), ("H", "O"))),
        ConnectivityTable({(0, 1), (1, 2)}),
    )
    assert thr.symbol() == "T"


def test_V():
    val = RGroup(
        group_set((("C", "H"), ("C", "H", "H", "H"), ("C", "H", "H", "H"))),
        ConnectivityTable({(0, 1), (1, 2)}),
    )
    assert val.symbol() == "V"


def test_W():
    trp = RGroup(
        group_set(
            (
                ("C",),
                ("C",),
                ("C",),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H", "H"),
                ("H", "N"),
            )
        ),
        ConnectivityTable({(0, 1)}),
    )
    assert trp.symbol() == "W"


def test_Y():
    tyr = RGroup(
        group_set(
            (
                ("C",),
                ("C",),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H"),
                ("C", "H", "H"),
                ("O", "H"),
            )
        ),
        ConnectivityTable({(0, 1)}),
    )
    assert tyr.symbol() == "Y"


def test_IL_invalid():
    with pytest.raises(KeyError) as msg:
        invalid_IL = RGroup(
            group_set(
                (
                    ("C", "H"),
                    ("C", "H", "H"),
                    ("C", "H", "H", "H"),
                    ("C", "H", "H", "H"),
                )
            ),
            ConnectivityTable({(0, 1), (2, 3), (1, 3)}),
        )
        assert not invalid_IL.symbol()
    assert (
        msg.exconly(True) == "KeyError: 'This R-group is I/L, but is neither "
        "I nor L.'"
    )
