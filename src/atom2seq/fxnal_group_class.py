from atom2seq.cluster_class import Cluster
from atom2seq.indexed_graph_class import IndexedGraph


class FxnalGroup(IndexedGraph):
    def __init__(
        self,
        clusters: set[Cluster],
        bonds: set[tuple[int]],
        parent: int = -1,
        idx: int = -1,
    ):
        super().__init__(clusters, bonds, parent, idx)
        self._symbol = ""

    def __repr__(self) -> str:
        return f"Group({self._vertices}, {self._edges}, {self._parent}, {self._idx})"  # noqa

    def _cleanup(self) -> None:
        super()._cleanup()
        self._update_symbol_determine_rep()

    def _print_symbol(self) -> None:
        print(f"in _print_symbol: {self._symbol=}")

    def _update_symbol_determine_rep(self) -> None:
        """Automatically assigns this functional group a symbol corresponding
        to one of five functional groups, or if it contains a single cluster,
        the symbol of that cluster. Also automatically sets the representative
        cluster if this functional group is one of the five supported
        functional groups."""
        # Creates a list of the symbols of all of the clusters in this
        # functional group.
        symbols = [cluster.get_symbol() for cluster in self.vertex_list()]
        print(f"{self=}, {symbols=}")
        # If there is one cluster in this functional group, they have the same
        # symbol.
        if len(self._vertices) == 1:
            print("len(self._vertices) == 1")
            self._symbol = self.vertex_list()[0].get_symbol()
            self._rep = 0
        # If there are two clusters in this functional group, the only FG we
        # care about is carbonyl, so we check for it.
        elif len(self._vertices) == 2:
            print("len(self._vertices) == 2")
            if ("C" in symbols) and ("O" in symbols):
                print('("C" in symbols) and ("O" in symbols)')
                carbon = self.vertex_list()[symbols.index("C")].get_idx()
                oxygen = self.vertex_list()[symbols.index("O")].get_idx()
                if self.check_edge((carbon, oxygen)):
                    print("self.check_edge((carbon, oxygen))")
                    self._symbol = "C=O"
                    print(f"{self.get_symbol()=}")
                    self._rep = carbon
        # If there are three clusters in this functional group, the only FGs we
        # care about are amide and carboxylic acid, so we check for those.
        elif len(self._vertices) == 3:
            # They both have a carbon-only cluster and an oxygen-only cluster
            # as they both contain a carbonyl.
            if ("C" in symbols) and ("O" in symbols):
                carbon = self.vertex_list()[symbols.index("C")].get_idx()
                oxygen = self.vertex_list()[symbols.index("O")].get_idx()
                if "NH2" in symbols:
                    amine = self.vertex_list()[symbols.index("NH2")].get_idx()
                    if self.check_edge((carbon, oxygen)) and self.check_edge(
                        (carbon, amine)
                    ):
                        self._symbol = "Amd"
                        self._rep = carbon
                elif "OH" in symbols:
                    hydroxyl = self.vertex_list()[symbols.index("OH")].get_idx()  # noqa
                    if self.check_edge((carbon, oxygen)) and self.check_edge(
                        (carbon, hydroxyl)
                    ):
                        self._symbol = "COOH"
                        self._rep = carbon
        # If there are six clusters in this functional group, the only FG we
        # care about is phenyl, so we check for it.
        elif len(self._vertices) == 6:
            if ("C" in symbols) and (symbols.count("CH") == 5):
                chs = []
                for cluster in self.vertex_list():
                    if cluster.get_symbol() == "CH":
                        chs.append(cluster.get_idx())
                carbon = self.vertex_list()[symbols.index("C")].get_idx()
                chch_bonds = 0
                chc_bonds = 0
                for ch1 in chs:
                    for ch2 in chs:
                        if self.check_edge((ch1, ch2)):
                            chch_bonds += 1
                    if self.check_edge((ch1, carbon)):
                        chc_bonds += 1
                # CH-CH bonds got double counted, so we divide by two.
                chch_bonds /= 2
                if (chch_bonds == 4) and (chc_bonds == 2):
                    self._symbol = "Ph"
                    self._rep = carbon
        # If there are seven clusters in this functional group, the only FG we
        # care about is phenol, so we check for it.
        elif len(self._vertices) == 7:
            if (
                ("OH" in symbols)
                and (symbols.count("C") == 2)
                and (symbols.count("CH") == 4)
            ):
                chs = []
                for cluster in self.vertex_list():
                    if cluster.get_symbol() == "CH":
                        chs.append(cluster.get_idx())
                carbons = []
                for cluster in self.vertex_list():
                    if cluster.get_symbol() == "C":
                        carbons.append(cluster.get_idx())
                hydroxyl = self.vertex_list()[symbols.index("OH")].get_idx()
                chch_bonds = 0
                chc_bonds = 0
                for ch1 in chs:
                    for ch2 in chs:
                        if self.check_edge((ch1, ch2)):
                            chch_bonds += 1
                    for c in carbons:
                        if self.check_edge((ch1, c)):
                            chc_bonds += 1
                coh_bonds = 0
                for c in carbons:
                    if self.check_edge((c, hydroxyl)):
                        coh_bonds += 1
                chch_bonds /= 2
                if (chch_bonds == 2) and (chc_bonds == 4) and (coh_bonds == 1):
                    self._symbol = "PhOH"
                    for carbon in carbons:
                        if self.check_edge((carbon, hydroxyl)):
                            self._rep = carbon
        else:
            # If this functional group has more than one cluster and isn't any
            # of our five functional groups, we assign it a standard name based
            # on the clusters inside of it. For example, an alkyne would be
            # (C)2 and a sulphonic acid would be (O)2(OH)(S).
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
        print(f"At the end of _update_symbol_determine_rep: {self._symbol=}")

    def get_symbol(self) -> str:
        """Returns the symbol of this functional group."""
        return self._symbol

    def get_atoms(self) -> set:
        atoms = set([])
        for cluster in self._vertices:
            atoms = atoms.union(cluster.get_vertices())
        return atoms

    def merge_clusters(self, idx1: int, idx2: int) -> None:
        if (idx1 in self.vertex_indices()) and (idx2 in self.vertex_indices()):
            if self.check_edge((idx1, idx2)):
                cluster1 = self.used_indices[idx1]
                cluster2 = self.used_indices[idx2]
                # Gets all atoms and bonds within these two clusters.
                atoms = cluster1.get_vertices().union(cluster2.get_vertices())
                bonds = cluster1.get_edges().union(cluster2.get_edges())
                # Adds the bond between the representative atoms of these
                # clusters.
                bonds.add((cluster1.get_rep(), cluster2.get_rep()))
                new_cluster = Cluster(atoms, bonds)
                # Finds all nodes clusters to these two clusters.
                bonded = self.get_adjacent(idx1).union(self.get_adjacent(idx2))
                bonded.remove(idx1)
                bonded.remove(idx2)
                # Deletes these two clusters.
                self.del_vertex(idx1)
                self.del_vertex(idx2)
                new_idx = new_cluster.get_idx()
                # Adds the new cluster and bonds it to all clusters that one of
                # the original two clusters was bonded to.
                self.add_vertex(new_cluster)
                print(f"in merge_clusters: {self._edges=}")
                print(self)
                for idx in bonded:
                    print(self.vertex_indices())
                    print(f"Bonding {idx} and {new_idx}")
                    self.add_edge((idx, new_idx))
            else:
                raise ValueError(
                    f"The objects at indices {idx1} and {idx2} are not bonded "
                    "and therefore cannot be merged."
                )
        else:
            raise ValueError(
                f"One or both of the objects at indices {idx1} or {idx2} are "
                "not children of this functional group."
            )

    def dist(self, other) -> float:
        self_rep = self.used_indices[self._rep]
        other_rep = self.used_indices[other.get_rep()]
        return self_rep.dist(other_rep)
