from atom2seq.cluster_class import Cluster
from atom2seq.indexed_graph_class import IndexedGraph


class FxnalGroup(IndexedGraph):
    def __repr__(self) -> str:
        return (
            f"Group({self.vertex_list()}, {sorted(list(self._edges))}, "
            f"{self._parent}, {self._idx})"
        )  # noqa

    def _cleanup(self) -> None:
        super()._cleanup()
        self._update_symbol_determine_rep()

    def _update_symbol_determine_rep(self) -> None:
        """Automatically assigns this functional group a symbol corresponding
        to one of five functional groups, or if it contains a single cluster,
        the symbol of that cluster. Also automatically sets the representative
        cluster if this functional group is one of the five supported
        functional groups."""
        # Creates a list of the symbols of all of the clusters in this
        # functional group.
        symbols = [cluster.get_symbol() for cluster in self.vertex_list()]
        # Creates a dictionary of what to call if there are a certain amount of
        # clusters. This way of doing a long if-elif-elif-... statement was
        # taken from https://stackoverflow.com/questions/17166074/ answered by
        # user Aya.
        call_this = {
            1: self._label_single,
            2: self._find_carbonyl,
            3: self._find_amide_or_cooh,
            6: self._find_phenyl,
            7: self._find_phenol,
        }
        if len(self._vertices) in call_this:
            call_this[len(self._vertices)](symbols)
        else:
            self._label_other(symbols)

    def _label_single(self, symbols):
        self._symbol = symbols[0]
        self._rep = 0

    def _find_carbonyl(self, symbols):
        if ("C" in symbols) and ("O" in symbols):
            carbon = self.vertex_list()[symbols.index("C")].get_idx()
            oxygen = self.vertex_list()[symbols.index("O")].get_idx()
            if self.check_edge((carbon, oxygen)):
                self._symbol = "C=O"
                self._rep = carbon

    def _find_amide_or_cooh(self, symbols):
        # Amide and COOH both have a carbon-only cluster and an oxygen-only
        # cluster, as they both contain a carbonyl.
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

    def _find_phenyl(self, symbols):
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
            # CH-CH bonds are being double counted
            chch_bonds /= 2
            if (chch_bonds == 4) and (chc_bonds == 2):
                self._symbol = "Ph"
                self._rep = carbon

    def _find_phenol(self, symbols):
        print("OH" in symbols, symbols.count("C"), symbols.count("CH"))
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
                    print(self.check_edge((ch1, c)))
                    if self.check_edge((ch1, c)):
                        chc_bonds += 1
            # CH-CH bonds are being double-counted
            chch_bonds /= 2
            coh_bonds = 0
            for c in carbons:
                if self.check_edge((c, hydroxyl)):
                    coh_bonds += 1
            print(chch_bonds, chc_bonds, coh_bonds)
            if (chch_bonds == 2) and (chc_bonds == 4) and (coh_bonds == 1):
                self._symbol = "PhOH"
                for carbon in carbons:
                    if not self.check_edge((carbon, hydroxyl)):
                        self._rep = carbon

    def _label_other(self, symbols):
        # If this functional group has more than one cluster and isn't any of
        # our five functional groups, we assign it a standard name based on the
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
        print(f"{self=}, {self._rep=}, {other=}, {other.get_rep()=}")
        print(
            "In Group.dist: Finding distance from "
            f"{self.used_indices[self_rep.get_rep()].coords} to "
            f"{self.used_indices[other_rep.get_rep()].coords}"
        )
        return self_rep.dist(other_rep)
