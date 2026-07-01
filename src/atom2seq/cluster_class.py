from atom2seq.indexed_graph_class import IndexedGraph


class Cluster(IndexedGraph):
    """A class to represent a cluster of atoms, typically a non-hydrogen
    surrounded by hydrogens."""

    def __repr__(self) -> str:
        return (
            f"Cluster({self.vertex_list()}, {sorted(list(self._edges))}, "
            f"{self._parent}, {self._idx})"
        )  # noqa

    def _tuple(self) -> tuple:
        """Returns a tuple representing this cluster. Used for checking
        equality and hashing."""
        return (
            self._symbol,
            tuple(sorted(tuple(self._vertices))),
            tuple(sorted(tuple(self._edges))),
        )  # noqa

    def _cleanup(self) -> None:
        super()._cleanup()
        self._update_representative()
        self._update_symbol()

    def _update_symbol(self) -> None:
        """Updates the symbol of this cluster by looking at the atoms inside of
        it."""
        # Initializes the output and a dictionary of number of given atoms with
        # certain symbols.
        symbol = ""
        num_syms = {}
        for atom in self._vertices:
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

    def _update_representative(self) -> None:
        """Sets the parent for all atoms within this cluster to this cluster's
        index. Also updates the representative atom for this cluster."""
        # Loops over all atoms, and if they are not a hydrogen, makes them the
        # cluster's new representative (for the intended use case, there will
        # only be one non-hydrogen atom.).
        for atom in self._vertices:
            if atom.symbol != "H":
                self._rep = atom.get_idx()

    def get_symbol(self) -> str:
        """Returns this cluster's symbol."""
        return self._symbol

    def dist(self, other) -> float:
        """Returns the distance from this cluster's representative atom to
        another cluster's representative atom."""
        self_rep = self.used_indices[self._rep]
        other_rep = self.used_indices[other.get_rep()]
        print(f"{self=}, {self._rep=}, {other=}, {other.get_rep()=}")
        print(
            f"In Cluster.dist: Finding distance from {self_rep.coords} to "
            f"{other_rep.coords}"
        )
        return self_rep.dist(other_rep)
