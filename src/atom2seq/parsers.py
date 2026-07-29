from time import time

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
    print("parsed")
    bond_mol(molecule)
    print("bonded")
    return molecule


def bond_mol(molecule: Mol) -> None:
    radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}
    for atom1 in molecule.atom_list():
        sym1 = atom1.symbol
        for atom2 in molecule.atom_list():
            sym2 = atom2.symbol
            max_dist = 1.15 * (radii[sym1] + radii[sym2])
            if (
                0.5
                <= molecule.dist(atom1.get_idx(), atom2.get_idx())
                <= max_dist  # noqa
            ):  # noqa
                molecule.get_bonds().add_pair(
                    (atom1.get_idx(), atom2.get_idx())
                )  # noqa


# def bond_mol(molecule: Mol) -> None:
#     radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}
#     data = KDTree(
#         np.array([list(atom.coords) for atom in molecule.atom_list()])
#     )  # noqa
#     bonds_by_idx = []
#     for atom in molecule.atom_list():
#         coords = list(atom.coords)
#         indices = data.query(coords, k=12)[1]
#         if isinstance(indices, np.ndarray):
#             bonds_by_idx.append(indices)
#         else:
#             bonds_by_idx.append([indices])
#     for i in range(len(bonds_by_idx)):
#         atom1 = False
#         for atom in molecule.atom_list():
#             if atom.get_idx() == i:
#                 atom1 = atom
#         bonds = bonds_by_idx[i]
#         for idx in bonds:
#             atom2 = False
#             for atom in molecule.atom_list():
#                 if atom.get_idx() == i:
#                     atom2 = atom
#             if idx < len(molecule.atom_list()):
#                 if (not molecule.get_bonds().check_pair((i, idx))) and (
#                     i != idx
#                 ):  # noqa
#                     if (
#                         0.5
#                         <= molecule.dist(i, idx)
#                         <= 1.15 * (radii[atom1.symbol] + radii[atom2.symbol])
#                     ):
#                         molecule.get_bonds().add_pair((i, idx))


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


def parse_nwc(filename: str) -> Mol:
    contents = file_base(filename)

    contents = [line.split() for line in contents]
    new_contents = []
    for line in contents:
        new_contents.append([line[1], *line[3:]])

    return parser_base(new_contents)
