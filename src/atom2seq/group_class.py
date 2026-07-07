import math

from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable


class Group:
    """A class representing a functional group. Supports checking equality."""

    idxs = {}

    def __init__(self, atoms: set[Atom], bonds: ConnectivityTable):
        self._atoms = atoms
        self._bonds = bonds
        self._idx = -1
        for i in range(len(self._atoms)):
            atom = sorted(list(self._atoms))[i]
            if atom.get_idx() == -1:
                atom.set_idx(i)

    def __repr__(self):
        return f"Group({self.atom_list()}, {self._bonds})"

    def __eq__(self, other):
        print(f"{(self._atoms == other.get_atoms())=}")
        print(f"{(self._bonds == other.get_bonds())=}")
        return (self._atoms == other.get_atoms()) and (
            self._bonds == other.get_bonds()
        )  # noqa

    def __hash__(self):
        return hash((self._atoms, self._bonds))

    def atom_list(self):
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

    def get_idx(self):
        """Returns the index of the group."""
        return self._idx

    def set_idx(self, new_idx):
        """Sets the index of the group to the given index. Also updates the
        internal dict of used indices."""
        old_idx = self._idx
        self._idx = new_idx
        if old_idx != -1:
            self.idxs.pop(old_idx)
        elif new_idx != -1:
            if new_idx not in self.idxs:
                self.idxs[new_idx] = self
            else:
                raise ValueError(
                    f"There is already an Atom at index {new_idx}."
                )  # noqa
