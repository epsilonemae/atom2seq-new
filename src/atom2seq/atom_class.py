class Atom:
    """A class representing an atom."""

    idxs = {}

    def __init__(self, symbol: str, coords: tuple[float], idx: int = -1):
        self.symbol = symbol
        self.coords = coords
        self._idx = idx
        if idx != -1:
            self.idxs[idx] = self

    def __repr__(self):
        return f"Atom('{self.symbol}', {self.coords})"

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
        """Sets the index of the atom to the given index. Also updates the
        internal dict of used indices."""
        old_idx = self._idx
        self._idx = new_idx
        if old_idx != -1:
            self.idxs.pop(old_idx)
        elif new_idx != -1:
            print(self.idxs)
            if new_idx in self.idxs:
                raise ValueError(
                    f"There is already an Atom at index {new_idx}."
                )  # noqa
            else:
                self.idxs[new_idx] = self
