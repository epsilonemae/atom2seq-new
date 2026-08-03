import math

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable


class Group:
    """A class representing a functional group. Supports checking equality."""

    def __init__(self, atoms: set[Atom], bonds: ConnectivityTable):
        self._atoms = atoms
        self._bonds = bonds
        self._idx = -1
        # Assigns indices to each of the atoms if and only if they didn't
        # already have an index
        for i in range(len(self._atoms)):
            atom = sorted(list(self._atoms))[i]
            if atom.get_idx() == -1:
                atom.set_idx(i)

    def __repr__(self):
        return f"Group({self.atom_list()}, {self._bonds}, {self._idx})"

    def __lt__(self, other):
        if len(self._atoms) == len(other.get_atoms()):
            if len(self._bonds) == len(other.get_bonds()):
                # If the lengths of the sets of atoms are the same and the
                # lengths of the ConnectivityTables of bonds are the same, then
                # the lt is based on the lt of the first atoms that is
                # different.
                for i in range(len(self._atoms)):
                    self_atom = self.atom_list()[i]
                    other_atom = other.atom_list()[i]
                    if self_atom.coords != other_atom.coords:
                        return self_atom < other_atom
            # If the lengths of the sets of atoms are the same, but the lengths
            # of the ConnectivityTables of bonds are different, the lt is based
            # on the lt of the lengths of the ConnectivityTables of bonds.
            else:
                return len(self._bonds) < len(other.get_bonds())
        # If the lengths of the sets of atoms are not the same, the lt is based
        # on the lt of those lengths
        else:
            return len(self._atoms) < len(other.get_atoms())

    def __eq__(self, other):
        return (self._atoms == other.get_atoms()) and (
            self._bonds == other.get_bonds()
        )  # noqa

    def __hash__(self):
        return hash((tuple(self.atom_list()), self._bonds))

    def atom_list(self) -> list:
        """Returns a sorted list of the atoms in the Group."""
        return sorted(list(self._atoms))

    def dist(self, n: int, m: int) -> float:
        """Calculates the Euclidean distance between the given atoms."""
        n_coords, m_coords = ((), ())
        for atom in self._atoms:
            if atom.get_idx() == n:
                n_coords = atom.coords
            elif atom.get_idx() == m:
                m_coords = atom.coords
        return math.sqrt(
            (n_coords[0] - m_coords[0]) ** 2
            + (n_coords[1] - m_coords[1]) ** 2
            + (n_coords[2] - m_coords[2]) ** 2
        )

    def get_bonds(self) -> ConnectivityTable:
        """Returns the ConnectivityTable of bonds."""
        return self._bonds

    def get_atoms(self) -> list[Atom]:
        """Returns the list of atoms."""
        return self._atoms

    def get_idx(self) -> int:
        """Returns the index of the group."""
        return self._idx

    def set_idx(self, new_idx: int) -> None:
        """Sets the index of the group to the given index."""
        self._idx = new_idx
