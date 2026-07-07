import math

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


class Mol:
    """A class representing a molecule. Supports checking equality."""

    def __init__(self, atoms: set[Atom], bonds: ConnectivityTable):
        self._atoms = atoms
        self._bonds = bonds
        for i in range(len(self._atoms)):
            self._atoms[i].set_idx(i)

    def __repr__(self):
        return f"Mol({self._atoms}, {self._bonds})"

    def __eq__(self, other):
        return (self._atoms == other.get_atoms()) and (
            self._bonds == other.get_bonds()
        )  # noqa

    def dist(self, n: int, m: int) -> float:
        """Calculates the Euclidean distance between the given atoms."""
        n_coords = self._atoms[n].coords
        m_coords = self._atoms[m].coords
        return math.sqrt(
            (n_coords[0] - m_coords[0]) ** 2
            + (n_coords[1] - m_coords[1]) ** 2
            + (n_coords[2] - m_coords[2]) ** 2
        )

    def get_bonds(self) -> ConnectivityTable:
        """Returns the list of bonds."""
        return self._bonds

    def get_atoms(self) -> set[Atom]:
        """Returns the list of atoms."""
        return self._atoms

    def group_atoms(self, idxs_to_group: list[int]) -> Group:
        """Returns a group object that contains the given atoms and the bonds
        between them."""
        atoms = set()
        bonds = ConnectivityTable([])
        for atom in self._atoms:
            if atom.get_idx() in idxs_to_group:
                atoms.add(atom)
        for bond in self._bonds:
            if (bond[0] in idxs_to_group) and (bond[1] in idxs_to_group):
                bonds.add_pair(bond)
        return Group(atoms, bonds)
