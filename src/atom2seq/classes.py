import math as m


class IndexedObject:
    """A class where each object within it must have a unique index. Supports
    checking equality, comparison, and hashing."""

    used_indices = {}

    def __init__(self, idx: int = -1):
        if idx == -1:
            idx = len(self.used_indices)
        if idx not in self.used_indices.keys():
            self._idx = idx
            self.used_indices[idx] = self
        else:
            raise ValueError(f"The index {idx} is already in use.")

    def __eq__(self, other) -> bool:
        # Indices are unique.
        if self._idx == other.get_idx():
            return True
        else:
            return self._idx_neq(other)

    def __lt__(self, other) -> bool:
        # Indices are unique.
        return self._idx < other.get_idx()

    def __hash__(self) -> int:
        return hash(self._tuple())

    def __repr__(self):
        return f"IndexedObject({self._idx})"

    def _tuple(self):
        return (self._idx,)

    def _idx_neq(self, other):
        """How to check if objects are equal if they have different indices."""
        return self._tuple() == other._tuple()

    def get_idx(self):
        return self._idx

    def set_idx(self, new_idx):
        if (new_idx not in self.used_indices.keys()) or (new_idx == self._idx):
            self.used_indices.pop(self._idx)
            self._idx = new_idx
            self.used_indices[new_idx] = self
        else:
            raise ValueError(f"The index {new_idx} is already in use.")


class Atom(IndexedObject):
    """A class to represent a single atom. Supports checking equality,
    comparison, and hashing."""

    def __init__(
        self, symbol: str, coords: tuple[float], parent: int = -1, idx: int = -1  # noqa
    ):  # noqa
        print(f"{symbol}, {coords}, {parent}")
        super().__init__(idx)
        self.symbol = symbol
        self.coords = coords
        self._parent = parent

    def __repr__(self):
        return (
            f"Atom('{self.symbol}', {self.coords}, {self._parent}, {self._idx})"  # noqa
        )

    def _tuple(self):
        return (self.symbol, self.coords)

    def get_idx(self) -> int:
        """Returns the index of the atom."""
        return self._idx

    def set_idx(self, new_idx: int) -> None:
        """Sets the index of the atom to the integer passed in."""
        self._idx = new_idx

    def get_parent(self) -> int:
        """Returns the index of this atom's parent."""
        return self._parent

    def set_parent(self, new_parent: int) -> None:
        """Sets the atom's parent to the object with the index passed in."""
        self._parent = new_parent

    def dist(self, other) -> float:
        """Returns the distance from this atom to another atom."""
        return m.sqrt(
            (self.coords[0] - other.coords[0]) ** 2
            + (self.coords[1] - other.coords[1]) ** 2
            + (self.coords[2] - other.coords[2]) ** 2
        )


class Node(IndexedObject):
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
        for atom in self._atoms:
            if atom.get_parent() != self._idx:
                atom.set_parent(self._idx)
            if atom.symbol != "H":
                self._rep = atom
        self._update_symbol()
        self._cleanup_bonds()

    def __hash__(self):
        return hash(
            (self._symbol, self._idx, tuple(self._atoms), tuple(self._bonds))
        )  # noqa

    def __repr__(self):
        return (
            f"Node({self._atoms}, {self._bonds}, {self._parent}, {self._idx})"  # noqa
        )

    def _tuple(self):
        return (
            self._symbol,
            sorted(tuple(self._atoms)),
            sorted(tuple(self._bonds)),
        )  # noqa

    def _update_symbol(self) -> None:
        symbol = ""
        num_syms = {}
        for atom in self._atoms:
            if atom.symbol in num_syms.keys():
                num_syms[atom.symbol] += 1
            else:
                num_syms[atom.symbol] = 1
        for sym in num_syms:
            if sym != "H":
                symbol += sym
                if num_syms[sym] != 1:
                    symbol += str(num_syms[sym])
        symbol += "H"
        if num_syms["H"] != 1:
            symbol += str(num_syms["H"])
        self._symbol = symbol

    def _verify_bonds(self, set_of_bonds):
        # Loops over all of the bonds
        for bond in set_of_bonds:
            print(f"{bond=}")
            found = [False, False]
            # Loops over each atom and checks if its index is one of the two in
            # the bond.
            for atom in self._atoms:
                print(f"{atom.get_idx()=}")
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
        new_bonds = set([])
        for bond in self._bonds:
            new_bonds.add((min(bond), max(bond)))

    def get_symbol(self) -> str:
        return self._symbol

    def get_atoms(self) -> set[Atom]:
        return self._atoms

    def set_atoms(self, new_atoms: set[Atom]) -> None:
        self._atoms = new_atoms
        self._bonds = {}

    def add_atom(self, symbol: str, coords: tuple[float]):
        self._atoms.add(Atom(symbol, coords, parent=self._idx))
        self._update_symbol()

    def del_atom(self, del_idx: int) -> None:
        for atom in self._atoms:
            if atom.get_idx() == del_idx:
                self._atoms.remove(atom)
                break
        self._update_symbol()

    def get_bonds(self) -> set[tuple[int]]:
        return self._bonds

    def set_bonds(self, new_bonds: set[tuple[int]]) -> None:
        if not self._verify_bonds(new_bonds):
            raise ValueError(
                "The bonds passed contain a bond involving an index not in this Node."  # noqa
            )
        else:
            self._bonds = new_bonds
            self._cleanup_bonds()

    def add_bond(self, idx1, idx2):
        if not self._verify_bonds({(idx1, idx2)}):
            raise ValueError(
                f"One of the indices {idx1} or {idx2} are not in this Node."
            )
        else:
            self._bonds.add((min(idx1, idx2), max(idx1, idx2)))

    def del_bond(self, idx1, idx2):
        if self.check_bond(idx1, idx2):
            self._bonds.remove((min(idx1, idx2), max(idx1, idx2)))

    def check_bond(self, idx1, idx2):
        for bond in self._bonds:
            if bond[0] == min(idx1, idx2) and bond[1] == max(idx1, idx2):
                return True
        return False

    def get_rep(self):
        return self._rep

    def set_rep(self, new_rep):
        if new_rep in {atom.get_idx() for atom in self._atoms}:
            self._rep = self.used_indices[new_rep]

    def get_parent(self):
        return self._parent

    def set_parent(self, new_parent):
        self._parent = new_parent

    def dist(self, other):
        return self.get_rep().dist(other.get_rep())
