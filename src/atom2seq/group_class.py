from atom2seq.cluster_class import Cluster
from atom2seq.indexed_object_class import IndexedObject


class Group(IndexedObject):
    def __init__(
        self,
        clusters: set[Cluster],
        bonds: set[tuple[int]],
        parent: int = -1,
        idx: int = -1,  # noqa
    ):  # noqa
        super().__init__(idx)
        self._clusters = clusters
        self._bonds = bonds
        self._parent = parent
        self._symbol = ""
        self._rep = False
        for cluster in self._clusters:
            cluster.set_parent(self._idx)
        self._update_symbol_determine_rep()
        self._cleanup_bonds()

    def _cleanup_bonds(self):
        """Sets each bond to be strictly increasing."""
        new_bonds = set([])
        for bond in self._bonds:
            new_bonds.add((min(bond), max(bond)))

    def _update_symbol_determine_rep(self):
        """Automatically assigns this group a symbol corresponding to one of
        five functional groups, or if it contains a single cluster, the symbol
        of that cluster. Also automatically sets the representative cluster if
        this group is one of the five supported functional groups."""
        # Creates a list of the symbols of all of the clusters in this group.
        symbols = [cluster.get_symbol() for cluster in self.cluster_list()]
        # If there is one cluster in this group, they have the same symbol.
        if len(self._clusters) == 1:
            self._symbol = self.cluster_list()[0].get_symbol()
            self._rep = self.cluster_list()[0]
        # If there are two clusters in this group, the only FG we care about is
        # carbonyl, so we check for it.
        elif len(self._clusters) == 2:
            if ("C" in symbols) and ("O" in symbols):
                carbon = self.cluster_list()[symbols.index("C")].get_idx()
                oxygen = self.cluster_list()[symbols.index("O")].get_idx()
                if self.check_bond(carbon, oxygen):
                    self._symbol = "C=O"
                    self._rep = self.cluster_list()[carbon]
        # If there are three clusters in this group, the only FGs we care about
        # are amide and carboxylic acid, so we check for those.
        elif len(self._clusters) == 3:
            # They both have a carbon-only cluster and an oxygen-only cluster
            # as they both contain a carbonyl.
            if ("C" in symbols) and ("O" in symbols):
                carbon = self.cluster_list()[symbols.index("C")].get_idx()
                oxygen = self.cluster_list()[symbols.index("O")].get_idx()
                if "NH2" in symbols:
                    amine = self.cluster_list()[symbols.index("NH2")].get_idx()
                    if self.check_bond(carbon, oxygen) and self.check_bond(
                        carbon, amine
                    ):
                        self._symbol = "Amd"
                        self._rep = self.cluster_list()[carbon]
                elif "OH" in symbols:
                    hydroxyl = self.cluster_list()[
                        symbols.index("OH")
                    ].get_idx()  # noqa
                    if self.check_bond(carbon, oxygen) and self.check_bond(
                        carbon, hydroxyl
                    ):
                        self._symbol = "COOH"
                        self._rep = self.cluster_list()[carbon]
        # If there are six clusters in this group, the only FG we care about is
        # phenyl, so we check for it.
        elif len(self._clusters) == 6:
            if ("C" in symbols) and (symbols.count("CH") == 5):
                chs = []
                for cluster in self.cluster_list():
                    if cluster.get_symbol() == "CH":
                        chs.append(cluster.get_idx())
                carbon = self.cluster_list()[symbols.index("C")].get_idx()
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
                    self._rep = self.cluster_list()[carbon]
        # If there are seven clusters in this group, the only FG we care about
        # is phenol, so we check for it.
        elif len(self._clusters) == 7:
            if (
                ("OH" in symbols)
                and (symbols.count("C") == 2)
                and (symbols.count("CH") == 4)
            ):
                chs = []
                for cluster in self.cluster_list():
                    if cluster.get_symbol() == "CH":
                        chs.append(cluster.get_idx())
                carbons = []
                for cluster in self.cluster_list():
                    if cluster.get_symbol() == "C":
                        carbons.append(cluster.get_idx())
                hydroxyl = self.cluster_list()[symbols.index("OH")].get_idx()
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
                    for carbon in carbons:
                        if self.check_bond(hydroxyl):
                            self._rep = self.cluster_list()[carbon]
        else:
            # If this group has more than one cluster and isn't any of our five
            # functional groups, we assign it a standard name based on the
            # clusters inside of it. For example, an alkyne would be (C)2 and a
            # sulphonic acid would be (O)2(OH)(S).
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

    def cluster_list(self):
        """Returns a sorted list copy of the set of clusters within this
        group."""
        return sorted(list(self._clusters))

    def get_symbol(self):
        """Returns the symbol of this group."""
        return self._symbol

    def check_bond(self, idx1, idx2):
        """Checks if two given atoms are bonded."""
        for bond in self._bonds:
            if bond[0] == min(idx1, idx2) and bond[1] == max(idx1, idx2):
                return True
        return False
