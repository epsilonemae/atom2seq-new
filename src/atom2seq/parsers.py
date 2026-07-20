import numpy as np
from scipy.spatial import KDTree

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.mol_class import Mol


def file_base(filename: str) -> list[str]:
    """Returns the stripped .readlines() of a passed file."""
    file = open(filename, "r")
    contents = file.readlines()
    file.close()
    contents = [line.strip() for line in contents]
    return contents


def parser_base(contents: list[list[str | int]]) -> Mol:
    """Takes in a list of lists of atomic symbols and their coordinates and
    returns a Mol containing all those atoms."""
    # Removes blank lines, changes any integers to be ints, and then returns a
    # Mol containing atoms made from each element of contents and no bonds.
    contents = [line for line in contents if line]
    new_contents = []
    for listy in contents:
        to_append = []
        for elt in listy:
            if elt.replace(".", "").replace("-", "").isdigit():
                to_append.append(float(elt))
            else:
                to_append.append(elt)
        new_contents.append(to_append)
    contents = new_contents
    atoms = set([Atom(listy[0], tuple(listy[1:])) for listy in contents])
    molecule = Mol(atoms, ConnectivityTable(set()))
    bond_mol(molecule)
    return molecule


def bond_mol(molecule: Mol) -> None:
    radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}
    max_bonds = {"H": 1, "O": 2, "N": 3, "C": 4, "S": 2}
    data = KDTree(
        np.array([list(atom.coords) for atom in molecule.get_atoms()])
    )  # noqa
    bonds_by_idx = []
    for atom in molecule.get_atoms():
        coords = list(atom.coords)
        indices = data.query(coords, k=max_bonds[atom.symbol])[1]
        if isinstance(indices, np.ndarray):
            bonds_by_idx.append([atom.get_idx(), indices])
        else:
            bonds_by_idx.append([atom.get_idx(), [indices]])
    for bond in bonds_by_idx:
        idx = bond[0]
        to_bond = bond[1]
        sym1 = molecule.get_atom(idx).symbol
        for idx_to_bond in to_bond:
            if idx_to_bond in molecule.idx_list():
                if idx != idx_to_bond:
                    sym2 = molecule.get_atom(idx_to_bond).symbol
                    max_dist = 1.1 * (radii[sym1] + radii[sym2])
                    if 0.5 <= molecule.dist(idx, idx_to_bond) <= max_dist:
                        molecule.get_bonds().add_pair((idx, idx_to_bond))


def parse_gjf(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .gjf format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    new_contents = ""
    for line in contents:
        new_contents += line
    contents = new_contents
    contents = contents.split("\\")
    contents = [line.split(",") for line in contents]

    return parser_base(contents)


def parse_xyz(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .xyz format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    # Checking if the first line is the number of atoms. If it is, remove
    # the first line and any blank lines that come after it.
    if contents[0][0].isdigit():
        contents.pop(0)
    contents = [line.split() for line in contents]

    return parser_base(contents)


def parse_pdb(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .pdb format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    contents = [line for line in contents if line[0:4] == "ATOM"]
    contents = [
        [line.split()[-1], *line.split()[-6:-3]]
        for line in contents
        if line.split()[-8] == "A"
    ]

    return parser_base(contents)


def parse_cif(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .cif format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    contents = [line for line in contents if line[0:4] == "ATOM"]
    contents = [
        [line.split()[2], *line.split()[10:13]]
        for line in contents
        if line.split()[6] == "A"
    ]

    return parser_base(contents)
