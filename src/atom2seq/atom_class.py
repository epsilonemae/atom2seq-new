class Atom:
    """A class representing an atom. Supports comparison, checking equality,
    and hashing."""

    def __init__(self, symbol: str, coords: tuple[float], idx: int = -1):
        self.symbol = symbol
        self.coords = coords
        self._idx = idx

    def __repr__(self):
        return f"Atom('{self.symbol}', {self.coords}, {self._idx})"

    def __eq__(self, other):
        # Returns True if and only if both the symbol and the coordinates are
        # the same
        return (self.symbol == other.symbol) and (self.coords == other.coords)

    def __lt__(self, other):
        if self.coords[0] == other.coords[0]:
            if self.coords[1] == other.coords[1]:
                return self.coords[2] < other.coords[2]
            else:
                return self.coords[1] < other.coords[1]
        else:
            return self.coords[0] < other.coords[0]

    def __hash__(self):
        return hash((self.symbol, self.coords))

    def get_idx(self):
        """Returns the index of the atom."""
        return self._idx

    def set_idx(self, new_idx):
        """Sets the index of the atom to the given index."""
        self._idx = new_idx
