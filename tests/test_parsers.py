import pytest

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.mol_class import Mol
from atom2seq.parsers import bond_mol  # noqa
from atom2seq.parsers import parse_cif, parse_gjf, parse_pdb, parse_xyz

prefix = __file__.removesuffix("test_parsers.py")


@pytest.fixture
def water_atoms():
    return {Atom("H", (1, 0, 0)), Atom("H", (0, 1, 0)), Atom("O", (0, 0, 0))}


def test_gjf_parser(water_atoms):
    out = parse_gjf(prefix + "assets/parser_tests/water.gjf")
    assert out.get_atoms() == water_atoms


def test_xyz_basic(water_atoms):
    out = parse_xyz(prefix + "assets/parser_tests/water.xyz")
    assert out.get_atoms() == water_atoms


def test_xyz_no_num(water_atoms):
    out = parse_xyz(prefix + "assets/parser_tests/water_no_number_of_atoms.xyz")  # noqa
    assert out.get_atoms() == water_atoms


def test_xyz_extra_lines(water_atoms):
    out = parse_xyz(prefix + "assets/parser_tests/water_extra_lines.xyz")
    assert out.get_atoms() == water_atoms


def test_xyz_no_lines(water_atoms):
    out = parse_xyz(prefix + "assets/parser_tests/water_no_lines.xyz")
    assert out.get_atoms() == water_atoms


def test_cif_parser(water_atoms):
    out = parse_cif(prefix + "assets/parser_tests/water.cif")
    assert out.get_atoms() == water_atoms


def test_cif_extra_info(water_atoms):
    out = parse_cif(prefix + "assets/parser_tests/water_extra_info.cif")
    assert out.get_atoms() == water_atoms


def test_pdb_parser(water_atoms):
    out = parse_pdb(prefix + "assets/parser_tests/water.pdb")
    assert out.get_atoms() == water_atoms


def test_pdb_extra_info(water_atoms):
    out = parse_pdb(prefix + "assets/parser_tests/water_extra_info.pdb")
    assert out.get_atoms() == water_atoms


radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}


@pytest.fixture
def pairs():
    return [
        ["H", "O"],
        ["H", "N"],
        ["H", "C"],
        ["C", "N"],
        ["C", "O"],
        ["S", "C"],
        ["S", "H"],
    ]


def out_of_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (1.2 * (radii[sym1] + radii[sym2]), 0, 0))
    return Mol({atom1, atom2}, ConnectivityTable(set()))


def at_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (1.1 * (radii[sym1] + radii[sym2]), 0, 0))
    return Mol({atom1, atom2}, ConnectivityTable(set()))


def in_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (radii[sym1] + radii[sym2], 0, 0))
    return Mol({atom1, atom2}, ConnectivityTable(set()))


def under_range(sym1, sym2):
    atom1 = Atom(sym1, (0, 0, 0))
    atom2 = Atom(sym2, (0.4, 0, 0))
    return Mol({atom1, atom2}, ConnectivityTable(set()))


def test_out_of_range(pairs):
    for pair in pairs:
        mol = out_of_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds().get_pairs() == set(), f"{pair}, {mol}"


def test_at_range(pairs):
    for pair in pairs:
        mol = at_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds().get_pairs() == {(0, 1)}


def test_in_range(pairs):
    for pair in pairs:
        mol = in_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds().get_pairs() == {(0, 1)}


def test_under_range(pairs):
    for pair in pairs:
        mol = under_range(*pair)
        bond_mol(mol)
        assert mol.get_bonds().get_pairs() == set()


@pytest.fixture
def glycine():
    return parse_xyz(prefix + "assets/glycine_optimized.xyz")


def test_glycine(glycine):
    glycine_bonds = {
        (5, 6),
        (4, 5),
        (4, 7),
        (0, 4),
        (0, 8),
        (0, 9),
        (1, 2),
        (1, 3),
        (0, 1),
    }
    glycine.get_bonds().get_pairs() == glycine_bonds
