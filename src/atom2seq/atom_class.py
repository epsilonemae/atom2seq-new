import math as m

from atom2seq.indexed_object_class import IndexedObject


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
