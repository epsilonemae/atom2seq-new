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
    """A class to represent a cluster of atoms, typically a non-hydrogen
    surrounded by hydrogens."""

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
        self._atom_parentage_and_representative()
        self._update_symbol()
        self._cleanup_bonds()

    def __repr__(self):
        return (
            f"Node({self._atoms}, {self._bonds}, {self._parent}, {self._idx})"  # noqa
        )

    def _tuple(self):
        """Returns a tuple representing this node. Used for checking equality
        and hashing."""
        return (
            self._symbol,
            tuple(sorted(tuple(self._atoms))),
            tuple(sorted(tuple(self._bonds))),
        )  # noqa

    def _update_symbol(self) -> None:
        """Updates the symbol of this node by looking at the atoms inside of
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

    def _verify_bonds(self, set_of_bonds):
        """Verifies that all of the bonds in a set are between atoms within
        this node."""
        # Loops over all of the bonds
        for bond in set_of_bonds:
            print(f"{bond=}")
            found = [False, False]
            # Loops over each atom and checks if its index is one of the two in
            # the bond.
            for atom in self._atoms:
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
        """Sets each bond to be strictly increasing."""
        new_bonds = set([])
        for bond in self._bonds:
            new_bonds.add((min(bond), max(bond)))

    def _atom_parentage_and_representative(self):
        """Sets the parent for all atoms within this node to this node's
        index. Also updates the representative atom for this node."""
        # Loops over all atoms, updates their parent, and if they are not a
        # hydrogen makes them this node's new representative (for the intended
        # use case, there will only be one non-hydrogen atom.).
        for atom in self._atoms:
            atom.set_parent(self._idx)
            if atom.symbol != "H":
                self._rep = atom

    def get_symbol(self) -> str:
        """Returns this node's symbol."""
        return self._symbol

    def get_atoms(self) -> set[Atom]:
        """Returns the list of atoms within this node."""
        return self._atoms

    def set_atoms(self, new_atoms: set[Atom]) -> None:
        """Sets the atoms within this node to the set passed in. In doing so,
        removes all bonds and updates the representative atom."""
        self._atoms = new_atoms
        self._bonds = {}
        self._atom_parentage_and_representative()

    def add_atom(self, symbol: str, coords: tuple[float]) -> None:
        """Adds a given atom to this node, updating the symbol as well."""
        self._atoms.add(Atom(symbol, coords, parent=self._idx))
        self._update_symbol()

    def del_atom(self, del_idx: int) -> None:
        """Deletes a given atom from this node, updating the symbol as well."""
        for atom in self._atoms:
            if atom.get_idx() == del_idx:
                self._atoms.remove(atom)
                break
        self._update_symbol()

    def get_bonds(self) -> set[tuple[int]]:
        """Returns the set of bonds within this node."""
        return self._bonds

    def set_bonds(self, new_bonds: set[tuple[int]]) -> None:
        """Sets the bonds within this node to the set passed in. Raises a
        ValueError if they are not all between atoms within this node."""
        if not self._verify_bonds(new_bonds):
            raise ValueError(
                "The bonds passed contain a bond involving an index not in this node."  # noqa
            )
        else:
            self._bonds = new_bonds
            self._cleanup_bonds()

    def add_bond(self, idx1, idx2):
        """Adds a bond to the current set of bonds. Raises a ValueError if it
        is not between two atoms within this node."""
        if not self._verify_bonds({(idx1, idx2)}):
            raise ValueError(
                f"One of the indices {idx1} or {idx2} are not in this node."
            )
        else:
            self._bonds.add((min(idx1, idx2), max(idx1, idx2)))

    def del_bond(self, idx1, idx2):
        """Deletes a given bond from the current set of bonds. Does nothing if
        that bond is not within this node."""
        if self.check_bond(idx1, idx2):
            self._bonds.remove((min(idx1, idx2), max(idx1, idx2)))

    def check_bond(self, idx1, idx2):
        """Checks if two given atoms are bonded."""
        for bond in self._bonds:
            if bond[0] == min(idx1, idx2) and bond[1] == max(idx1, idx2):
                return True
        return False

    def get_rep(self):
        """Returns the representative atom of this node."""
        return self._rep

    def set_rep(self, new_rep):
        """Sets the representative atom of this node to the atom with the given
        index."""
        if new_rep in {atom.get_idx() for atom in self._atoms}:
            self._rep = self.used_indices[new_rep]
        else:
            raise ValueError(
                f"The object at index {new_rep} is not an atom within this node."  # noqa
            )

    def get_parent(self):
        """Returns the index of this node's parent."""
        return self._parent

    def set_parent(self, new_parent):
        """Sets this node's parent to the object with the given index."""
        self._parent = new_parent

    def dist(self, other):
        """Returns the distance from this node's representative atom to another
        node's representative atom."""
        return self.get_rep().dist(other.get_rep())


class Group(IndexedObject):
    def __init__(
        self,
        nodes: set[Node],
        bonds: set[tuple[int]],
        parent: int = -1,
        idx: int = -1,  # noqa
    ):  # noqa
        super().__init__(idx)
        self._nodes = nodes
        self._bonds = bonds
        self._parent = parent
        for node in self._nodes:
            node.set_parent(self._idx)
        self._update_symbol()
        self._cleanup_bonds()

    def _cleanup_bonds(self):
        """Sets each bond to be strictly increasing."""
        new_bonds = set([])
        for bond in self._bonds:
            new_bonds.add((min(bond), max(bond)))

    def _update_symbol(self):
        """Automatically assigns this group a symbol corresponding to one of
        five functional groups, or if it contains a single node, the symbol of
        that node."""
        # Creates a list of the symbols of all of the nodes in this group.
        symbols = [node.get_symbol() for node in self.node_list()]
        # If there is one node in the group, they have the same symbol.
        if len(self._nodes) == 1:
            self._symbol = self.node_list()[0].get_symbol()
        # If there are two nodes in the group, the only
        elif len(self._nodes) == 2:
            if ("C" in symbols) and ("O" in symbols):
                carbon = self.node_list()[symbols.index("C")].get_idx()
                oxygen = self.node_list()[symbols.index("O")].get_idx()
                if self.check_bond(carbon, oxygen):
                    self._symbol = "C=O"
        elif len(self._nodes) == 3:
            if ("C" in symbols) and ("O" in symbols):
                carbon = self.node_list()[symbols.index("C")].get_idx()
                oxygen = self.node_list()[symbols.index("O")].get_idx()
                if "NH2" in symbols:
                    amine = self.node_list()[symbols.index("NH2")].get_idx()
                    if self.check_bond(carbon, oxygen) and self.check_bond(
                        carbon, amine
                    ):
                        self._symbol = "Amd"
                elif "OH" in symbols:
                    hydroxyl = self.node_list()[symbols.index("OH")].get_idx()
                    if self.check_bond(carbon, oxygen) and self.check_bond(
                        carbon, hydroxyl
                    ):
                        self._symbol = "COOH"
        elif len(self._nodes) == 6:
            if ("C" in symbols) and (symbols.count("CH") == 5):
                chs = []
                for node in self.node_list():
                    if node.get_symbol() == "CH":
                        chs.append(node.get_idx())
                carbon = self.node_list()[symbols.index("C")].get_idx()
                chch_bonds = 0
                chc_bonds = 0
                for ch1 in chs:
                    for ch2 in chs:
                        if self.check_bond(ch1, ch2):
                            chch_bonds += 1
                    if self.check_bond(ch1, carbon):
                        chc_bonds += 1
                # CH-CH bonds got double counted, so we divide by two.
                chch_bonds /= 2
                if (chch_bonds == 4) and (chc_bonds == 2):
                    self._symbol = "Ph"
        elif len(self._nodes) == 7:
            if (
                ("OH" in symbols)
                and (symbols.count("C") == 2)
                and (symbols.count("CH") == 4)
            ):
                chs = []
                for node in self.node_list():
                    if node.get_symbol() == "CH":
                        chs.append(node.get_idx())
                carbons = []
                for node in self.node_list():
                    if node.get_symbol() == "C":
                        carbons.append(node.get_idx())
                hydroxyl = self.node_list()[symbols.index("OH")].get_idx()
                chch_bonds = 0
                chc_bonds = 0
                for ch1 in chs:
                    for ch2 in chs:
                        if self.check_bond(ch1, ch2):
                            chch_bonds += 1
                    for c in carbons:
                        if self.check_bond(ch1, c):
                            chc_bonds += 1
                coh_bonds = 0
                for c in carbons:
                    if self.check_bond(c, hydroxyl):
                        coh_bonds += 1
                chch_bonds /= 2
                if (chch_bonds == 2) and (chc_bonds == 4) and (coh_bonds == 1):
                    self._symbol = "PhOH"
        else:
            out = ""
            num_syms = {}
            for sym in symbols:
                if sym in num_syms.keys():
                    num_syms[sym] = 1
                else:
                    num_syms[sym] += 1
            for sym in num_syms.keys():
                out += f"({sym}){num_syms[sym]}"
            self._symbol = out

    def node_list(self):
        return sorted(list(self._nodes))

    def get_symbol(self):
        return self._symbol

    def check_bond(self, idx1, idx2):
        """Checks if two given atoms are bonded."""
        for bond in self._bonds:
            if bond[0] == min(idx1, idx2) and bond[1] == max(idx1, idx2):
                return True
        return False
