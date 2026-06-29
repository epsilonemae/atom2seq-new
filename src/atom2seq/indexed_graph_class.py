from atom2seq.indexed_object_class import IndexedObject


class IndexedGraph(IndexedObject):
    """A class to represent a graph. This graph has an index and its vertices
    are indexed objects as well."""

    def __init__(
        self,
        vertices: set[IndexedObject],
        edges: set[tuple[int]],
        parent: int = -1,
        idx: int = -1,
    ):
        super().__init__(parent, idx)
        self._vertices = vertices
        self._edges = edges
        self._rep = -1
        self._cleanup()

    def __repr__(self):
        return (
            f"IndexedGraph({self._vertices}, {self._edges}, {self._parent}, "
            f"{self._idx})"
        )

    def _cleanup(self):
        for vertex in self._vertices:
            vertex.set_parent(self._idx)
        new_edges = set([])
        for edge in self._edges:
            if (edge[0] not in self.vertex_indices()) or (
                edge[1] not in self.vertex_indices()
            ):
                raise ValueError(
                    f"The edge {edge} contains the index of an object that is"
                    "not a child of this graph."
                )
            new_edges.add((min(edge), max(edge)))
        self._edges = new_edges

    def _tuple(self):
        return (
            tuple(sorted(tuple(self._vertices))),
            tuple(sorted(tuple(self._edges))),
        )  # noqa

    def vertex_indices(self):
        return {vertex.get_idx() for vertex in self._vertices}

    def vertex_list(self):
        return sorted(list(self._vertices))

    def get_vertices(self):
        return self._vertices

    def set_vertices(self, new_vertices: set[IndexedObject]):
        """Sets this graph's vertex set to the given set. Also removes all
        edges."""
        self._vertices = new_vertices
        self._edges = set([])
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
            for edge in self._edges:
                if (edge[0] == idx_to_del) or (edge[1] == idx_to_del):
                    self._edges.remove(edge)
            self._cleanup()

    def get_edges(self):
        return self._edges

    def set_edges(self, new_edges: set[tuple[int]]):
        self._edges = new_edges
        self._cleanup()

    def add_edge(self, new_edge: tuple[int]):
        self._edges.add(new_edge)
        self._cleanup()

    def del_edge(self, edge_to_del):
        if edge_to_del in self._edges:
            self._edges.remove(edge_to_del)
            self._cleanup()

    def check_edge(self, edge_to_check):
        return edge_to_check in self._edges

    def get_rep(self):
        return self._rep

    def set_rep(self, idx):
        if idx in self.vertex_indices():
            self._rep = idx
        else:
            raise ValueError(
                f"The object at index {idx} is not a child of this graph."
            )  # noqa

    def get_adjacent(self, idx):
        out = {}
        if idx in self.vertex_indices():
            for edge in self._edges:
                if edge[0] == idx:
                    out.add(edge[1])
                elif edge[1] == idx:
                    out.add(edge[0])
        return out
