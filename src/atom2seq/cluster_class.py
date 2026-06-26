from atom2seq.atom_class import Atom
from atom2seq.indexed_object_class import IndexedObject


class Cluster(IndexedObject):
    """A class to represent a cluster of atoms, typically a non-hydrogen
    surrounded by hydrogens."""

    def __init__(
        self,
        atoms: set[Atom],
        bonds: set[tuple[int]],
        parent: int = -1,
        idx: int = -1,  # noqa
    ):  # noqa
        super().__init__(idx)
        self._atoms = atoms
        self._bonds = bonds
        self._parent = parent
        self._atom_parentage_and_representative()
        self._update_symbol()
        self._cleanup_bonds()

    def __repr__(self):
        return f"Cluster({self._atoms}, {self._bonds}, {self._parent}, {self._idx})"  # noqa

    def _tuple(self):
        """Returns a tuple representing this cluster. Used for checking
        equality and hashing."""
        return (
            self._symbol,
            tuple(sorted(tuple(self._atoms))),
            tuple(sorted(tuple(self._bonds))),
        )  # noqa

    def _update_symbol(self) -> None:
        """Updates the symbol of this cluster by looking at the atoms inside of
        it."""
        # Initializes the output and a dictionary of number of given atoms with
        # certain symbols.
        symbol = ""
        num_syms = {}
        for atom in self._atoms:
            # Adds one to the current value for this atom's symbol (if it
            # doesn't have one, sets it to one)
            if atom.symbol in num_syms.keys():
                num_syms[atom.symbol] += 1
            else:
                num_syms[atom.symbol] = 1
        # Loops over all symbols in the dictionary and puts all non-hydrogen
        # atoms first (e.g. C2N3...)
        for sym in num_syms:
            if sym != "H":
                symbol += sym
                if num_syms[sym] != 1:
                    symbol += str(num_syms[sym])
        # Now adds the hydrogens (e.g. C2N3H5)
        if "H" in num_syms.keys():
            symbol += "H"
            if num_syms["H"] != 1:
                symbol += str(num_syms["H"])
        self._symbol = symbol

    def _verify_bonds(self, set_of_bonds):
        """Verifies that all of the bonds in a set are between atoms within
        this cluster."""
        # Loops over all of the bonds
        for bond in set_of_bonds:
            print(f"{bond=}")
            found = [False, False]
            # Loops over each atom and checks if its index is one of the two in
            # the bond.
            for atom in self._atoms:
                if atom.get_idx() == bond[0]:
                    found[0] = True
                elif atom.get_idx() == bond[1]:
                    found[1] = True
                if found[0] and found[1]:
                    break
            if (not found[0]) or (not found[1]):
                return False
        return True

    def _cleanup_bonds(self):
        """Sets each bond to be strictly increasing."""
        new_bonds = set([])
        for bond in self._bonds:
            new_bonds.add((min(bond), max(bond)))

    def _atom_parentage_and_representative(self):
        """Sets the parent for all atoms within this cluster to this cluster's
        index. Also updates the representative atom for this cluster."""
        # Loops over all atoms, updates their parent, and if they are not a
        # hydrogen makes them thcluster's new representative (for the intended
        # use case, there will only be one non-hydrogen atom.).
        for atom in self._atoms:
            atom.set_parent(self._idx)
            if atom.symbol != "H":
                self._rep = atom

    def get_symbol(self) -> str:
        """Returns this cluster's symbol."""
        return self._symbol

    def get_atoms(self) -> set[Atom]:
        """Returns the list of atoms within this cluster."""
        return self._atoms

    def set_atoms(self, new_atoms: set[Atom]) -> None:
        """Sets the atoms within this cluster to the set passed in. In doing
        so, removes all bonds and updates the representative atom."""
        self._atoms = new_atoms
        self._bonds = {}
        self._atom_parentage_and_representative()

    def add_atom(self, symbol: str, coords: tuple[float]) -> None:
        """Adds a given atom to this cluster, updating the symbol as well."""
        self._atoms.add(Atom(symbol, coords, parent=self._idx))
        self._update_symbol()

    def del_atom(self, del_idx: int) -> None:
        """Deletes a given atom from this cluster, updating the symbol as
        well."""
        for atom in self._atoms:
            if atom.get_idx() == del_idx:
                self._atoms.remove(atom)
                break
        self._update_symbol()

    def get_bonds(self) -> set[tuple[int]]:
        """Returns the set of bonds within this cluster."""
        return self._bonds

    def set_bonds(self, new_bonds: set[tuple[int]]) -> None:
        """Sets the bonds within this cluster to the set passed in. Raises a
        ValueError if they are not all between atoms within this cluster."""
        if not self._verify_bonds(new_bonds):
            raise ValueError(
                "The bonds passed contain a bond involving an index not in this cluster."  # noqa
            )
        else:
            self._bonds = new_bonds
            self._cleanup_bonds()

    def add_bond(self, idx1, idx2):
        """Adds a bond to the current set of bonds. Raises a ValueError if it
        is not between two atoms within this cluster."""
        if not self._verify_bonds({(idx1, idx2)}):
            raise ValueError(
                f"One of the indices {idx1} or {idx2} are not in this cluster."
            )
        else:
            self._bonds.add((min(idx1, idx2), max(idx1, idx2)))

    def del_bond(self, idx1, idx2):
        """Deletes a given bond from the current set of bonds. Does nothing if
        that bond is not within this cluster."""
        if self.check_bond(idx1, idx2):
            self._bonds.remove((min(idx1, idx2), max(idx1, idx2)))

    def check_bond(self, idx1, idx2):
        """Checks if two given atoms are bonded."""
        for bond in self._bonds:
            if bond[0] == min(idx1, idx2) and bond[1] == max(idx1, idx2):
                return True
        return False

    def get_rep(self):
        """Returns the representative atom of this cluster."""
        return self._rep

    def set_rep(self, new_rep):
        """Sets the representative atom of this cluster to the atom with the
        given index."""
        if new_rep in {atom.get_idx() for atom in self._atoms}:
            self._rep = self.used_indices[new_rep]
        else:
            raise ValueError(
                f"The object at index {new_rep} is not an atom within this cluster."  # noqa
            )

    def get_parent(self):
        """Returns the index of this cluster's parent."""
        return self._parent

    def set_parent(self, new_parent):
        """Sets this cluster's parent to the object with the given index."""
        self._parent = new_parent

    def dist(self, other):
        """Returns the distance from this cluster's representative atom to
        another cluster's representative atom."""
        return self.get_rep().dist(other.get_rep())
