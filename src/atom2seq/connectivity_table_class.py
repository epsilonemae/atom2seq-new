class ConnectivityTable:
    """A class to represent the edges of a graph. The vertices are not stored
    within this class."""

    def __init__(self, pairs: set[tuple[int]]):
        self._pairs = pairs

    def __eq__(self, other):
        return self._pairs == other.get_pairs()

    def __repr__(self):
        return f"ConnectivityTable({self._pairs})"

    def __len__(self):
        return len(self._pairs)

    def __hash__(self):
        return hash(tuple(self._pairs))

    def get_pairs(self):
        """Returns the list of pairs within this table."""
        return self._pairs

    def add_pair(self, new_pair: tuple[int]):
        """Adds a given pair to this table."""
        self._pairs.add(new_pair)

    def del_pair(self, pair_to_del: tuple[int]):
        """Deletes a given pair from this table."""
        self._pairs.remove(pair_to_del)

    def check_pair(self, pair_to_check: tuple[int]):
        """Checks if a given pair is in this table."""
        return pair_to_check in self._pairs

    def get_paired(self, index):
        """Returns a list of all indices paired with a given index."""
        out = set()
        for pair in self._pairs:
            if pair[0] == index:
                out.add(pair[1])
            elif pair[1] == index:
                out.add(pair[0])
        return out
