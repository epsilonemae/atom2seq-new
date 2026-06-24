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
        return hash((self._idx, *self._hash_tuple()))

    def __repr__(self):
        return f"IndexedObject({self._idx})"

    def _hash_tuple(self):
        return (self._idx,)

    def _idx_neq(self, other):
        """How to check if objects are equal if they have different indices."""
        return self._hash_tuple() == other._hash_tuple()

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
        self, symbol: str, coords: tuple[float], parent: int, idx: int = -1
    ):  # noqa
        super().__init__(idx)
        self.symbol = symbol
        self.coords = coords
        self._parent = parent

    def __repr__(self):
        return (
            f"Atom('{self.symbol}', {self.coords}, {self._idx}, {self._parent})"  # noqa
        )

    def _hash_tuple(self):
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
        """Sets the atom's parent to the parent with the index passed in."""
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
        self, atoms: set[Atom], bonds: set[tuple], idx: int, parent: int
    ):  # noqa
        super().__init__(idx)
        self._atoms = atoms
        self._bonds = bonds
        self._parent = parent
        self._update_symbol()

    def __hash__(self):
        return hash(
            (self._symbol, self._idx, tuple(self._atoms), tuple(self._bonds))
        )  # noqa

    def _hash_tuple(self):
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
                symbol += sym + num_syms[sym]
        symbol += "H" + num_syms["H"]
        self._symbol = symbol

    def get_symbol(self) -> str:
        return self._symbol

    def get_atoms(self) -> set[Atom]:
        return self._atoms

    def set_atoms(self, new_atoms: set[Atom]) -> None:
        self._atoms = new_atoms

    def add_atom(self, symbol: str, coords: tuple[float]):
        self._atoms.add(Atom(symbol, coords, len(self._atoms), self._idx))
        self._update_symbol()

    def del_atom(self, del_idx: int) -> None:
        for atom in self._atoms:
            if atom.get_idx() == del_idx:
                self._atoms.remove(atom)
                break
        self._update_symbol()
