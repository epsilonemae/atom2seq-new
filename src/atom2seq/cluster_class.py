from atom2seq.indexed_graph_class import IndexedGraph


class Cluster(IndexedGraph):
    """A class to represent a cluster of atoms, typically a non-hydrogen
    surrounded by hydrogens."""

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

    def _cleanup(self):
        super()._cleanup()
        self._atom_parentage_and_representative()
        self._update_symbol()

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

    def dist(self, other):
        """Returns the distance from this cluster's representative atom to
        another cluster's representative atom."""
        self_rep = self.used_indices[self._rep]
        other_rep = self.used_indices[other.get_rep()]
        return self_rep.dist(other_rep)
