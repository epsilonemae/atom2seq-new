import math as m

from atom2seq.indexed_object_class import IndexedObject


class Atom(IndexedObject):
    """A class to represent a single atom. Supports checking equality,
    comparison, and hashing."""

    def __init__(
        self, symbol: str, coords: tuple[float], parent: int = -1, idx: int = -1  # noqa
    ):  # noqa
        super().__init__(idx)
        self.symbol = symbol
        self.coords = coords

    def __repr__(self) -> str:
        return (
            f"Atom('{self.symbol}', {self.coords}, {self._parent}, {self._idx})"  # noqa
        )

    def _tuple(self) -> tuple:
        return (self.symbol, self.coords)

    def dist(self, other) -> float:
        """Returns the distance from this atom to another atom."""
        print(f"{self=}, {self.coords=}, {other}, {other.coords}")
        print(
            f"In Atom.dist: Finding the distance from {self.coords} to "
            f"{other.coords}"
        )
        return m.sqrt(
            (self.coords[0] - other.coords[0]) ** 2
            + (self.coords[1] - other.coords[1]) ** 2
            + (self.coords[2] - other.coords[2]) ** 2
        )
