import numpy as np
from scipy.spatial import KDTree

from atom2seq.classes import Mol

radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}

max_bonds = {"H": 1, "O": 2, "N": 3, "C": 4, "S": 2}


def bond_mol(molecule: Mol) -> None:
    data = KDTree(
        np.array([list(atom.coords) for atom in molecule.get_atoms()])
    )  # noqa
    bonds_by_idx = []
    for atom in molecule.get_atoms():
        # print(atom)
        coords = list(atom.coords)
        indices = data.query(coords, k=max_bonds[atom.symbol])[1]
        # print(type(indices), indices)
        if isinstance(indices, np.ndarray):
            bonds_by_idx.append(indices)
        else:
            bonds_by_idx.append([indices])
    for i in range(len(bonds_by_idx)):
        bonds = bonds_by_idx[i]
        for idx in bonds:
            if idx < len(molecule.get_atoms()):
                if (not molecule.is_bond(i, idx)) and (i != idx):
                    if (
                        0.5
                        <= molecule.dist(i, idx)
                        <= 1.1
                        * (
                            radii[molecule.get_atoms()[i].symbol]
                            + radii[molecule.get_atoms()[idx].symbol]
                        )
                    ):
                        molecule.add_bond(i, idx)
