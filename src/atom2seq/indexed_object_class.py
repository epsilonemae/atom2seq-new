class IndexedObject:
    """A class where each object within it must have a unique index. Supports
    checking equality, comparison, and hashing."""

    used_indices = {}

    def __init__(self, parent: int = -1, idx: int = -1):
        self._parent = parent
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

    def _reset_used_indices(self):
        to_pop = []
        for elt in self.used_indices:
            to_pop.append(elt)
        for elt in to_pop:
            self.used_indices.pop(elt)
        print(f"in _reset_used_indices: {self.used_indices=}")

    def get_idx(self):
        return self._idx

    def set_idx(self, new_idx):
        if (new_idx not in self.used_indices.keys()) or (new_idx == self._idx):
            self.used_indices.pop(self._idx)
            self._idx = new_idx
            self.used_indices[new_idx] = self
        else:
            raise ValueError(f"The index {new_idx} is already in use.")

    def get_parent(self):
        return self.used_indices[self._parent]

    def set_parent(self, new_parent):
        if new_parent in self.used_indices.keys():
            self._parent = new_parent
        else:
            raise ValueError(
                f"There is no object with index {new_parent}, so it cannot "
                f"be assigned as the object with index {self._idx}'s parent."  # noqa
            )
