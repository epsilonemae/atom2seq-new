from atom2seq.indexed_object_class import IndexedObject


class IndexedGraph(IndexedObject):
    """A class to represent a graph. This graph has an index and its vertices
    are indexed objects as well."""

    def __init__(
        self,
        vertices: set[IndexedObject],
        bonds: set[tuple[int]],
        parent: int = -1,
        idx: int = -1,
    ):
        super().__init__(parent, idx)
        self._vertices = vertices
        self._bonds = bonds
        self._cleanup()

    def __repr__(self):
        return (
            f"IndexedGraph({self._vertices}, {self._bonds}, {self._parent}, "
            f"{self._idx})"
        )

    def _cleanup(self):
        for vertex in self._vertices:
            vertex.set_parent(self._idx)
        new_bonds = {}
        for bond in self._bonds:
            if (bond[0] not in self.vertex_indices()) or (
                bond[1] not in self.vertex_indices()
            ):
                raise ValueError(
                    f"The bond {bond} contains the index of an object that is"
                    "not a child of this graph."
                )
            new_bonds.add((min(bond), max(bond)))
        self._bonds = new_bonds

    def _tuple(self):
        return (
            tuple(sorted(tuple(self._vertices))),
            tuple(sorted(tuple(self._bonds))),
        )  # noqa

    def vertex_indices(self):
        return {vertex.get_idx() for vertex in self._vertices}

    def vertex_list(self):
        return sorted(list(self._vertices))

    def get_vertices(self):
        return self._vertices

    def set_vertices(self, new_vertices: set[IndexedObject]):
        """Sets this graph's vertex set to the given set. Also removes all
        bonds."""
        self._vertices = new_vertices
        self._bonds = set([])
        self._cleanup()

    def add_vertex(self, new_vertex: IndexedObject):
        self._vertices.add(new_vertex)
        self._cleanup()

    def del_vertex(self, idx_to_del: int):
        if idx_to_del in self.vertex_indices():
            to_del = ""
            for vertex in self._vertices:
                if vertex.get_idx() == idx_to_del:
                    to_del = vertex
                    break
            self._vertices.remove(to_del)
            for bond in self._bonds:
                if (bond[0] == idx_to_del) or (bond[1] == idx_to_del):
                    self._bonds.remove(bond)
            self._cleanup()

    def get_bonds(self):
        return self._bonds

    def set_bonds(self, new_bonds: set[tuple[int]]):
        self._bonds = new_bonds
        self._cleanup()

    def add_bond(self, new_bond: tuple[int]):
        self._bonds.add(new_bond)
        self._cleanup()

    def del_bond(self, bond_to_del):
        if bond_to_del in self._bonds:
            self._bonds.remove(bond_to_del)
            self._cleanup()

    def check_bond(self, bond_to_check):
        return bond_to_check in self._bonds
