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
            self.atom_list()[i].set_idx(i)

    def __repr__(self):
        return f"Mol({self.atom_list()}, {self._bonds})"

    def __eq__(self, other):
        return (self._atoms == other.get_atoms()) and (
            self._bonds == other.get_bonds()
        )  # noqa

    def atom_list(self):
        return sorted(list(self._atoms))

    def idx_list(self):
        return [atom.get_idx() for atom in self.atom_list()]

    def get_atom(self, idx):
        """Returns the atom at the specified index."""
        return self.atom_list()[self.idx_list().index(idx)]

    def dist(self, n: int, m: int) -> float:
        """Calculates the Euclidean distance between the given atoms."""
        n_coords, m_coords = ((), ())
        n_coords = self.get_atom(n).coords
        m_coords = self.get_atom(m).coords
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
        new_atoms = set()
        bonds = ConnectivityTable(set())
        for atom in self._atoms:
            if atom.get_idx() in idxs_to_group:
                new_atoms.add(atom)
        for bond in self._bonds._pairs:
            if (bond[0] in idxs_to_group) and (bond[1] in idxs_to_group):
                bonds.add_pair(bond)
        return Group(new_atoms, bonds)
