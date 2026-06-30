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

    def __repr__(self) -> str:
        return (
            f"IndexedGraph({self._vertices}, {self._edges}, {self._parent}, "
            f"{self._idx})"
        )

    def _cleanup(self) -> None:
        """Performs various cleanup tasks, such as setting the parent for all
        of the vertices in this graph to this graph, making sure all edges are
        within this graph, and making edges always have the lower index
        first."""
        for vertex in self._vertices:
            vertex.set_parent(self._idx)
        new_edges = set([])
        print(f"in super()._cleanup: {self._edges=}, {self=}")
        for edge in self._edges:
            if (edge[0] not in self.vertex_indices()) or (
                edge[1] not in self.vertex_indices()
            ):
                raise ValueError(
                    f"The edge {edge} contains the index of an object that is "
                    "not a child of this graph."
                )
            new_edges.add((min(edge), max(edge)))
        self._edges = new_edges

    def _tuple(self) -> tuple:
        return (
            tuple(sorted(tuple(self._vertices))),
            tuple(sorted(self._edges)),
        )  # noqa

    def vertex_indices(self) -> set:
        """Returns a set containing the indices of all vertices within this
        graph."""
        return {vertex.get_idx() for vertex in self._vertices}

    def vertex_list(self) -> list:
        """Returns a sorted list of the vertices within this graph."""
        return sorted(list(self._vertices))

    def get_vertices(self) -> set:
        """Returns the set of vertices within this graph."""
        return self._vertices

    def set_vertices(self, new_vertices: set[IndexedObject]) -> None:
        """Sets this graph's vertex set to the given set. Also removes all
        edges."""
        self._vertices = new_vertices
        self._edges = set([])
        self._cleanup()

    def add_vertex(self, new_vertex: IndexedObject) -> None:
        """Adds a given vertex to this graph."""
        self._vertices.add(new_vertex)
        self._cleanup()

    def del_vertex(self, idx_to_del: int) -> None:
        """Deletes the vertex at the given index from this graph. Does not
        delete the vertex object itself."""
        if idx_to_del in self.vertex_indices():
            to_del = ""
            for vertex in self._vertices:
                if vertex.get_idx() == idx_to_del:
                    to_del = vertex
                    break
            self._vertices.remove(to_del)
            edges_to_del = set([])
            for edge in self._edges:
                if (edge[0] == idx_to_del) or (edge[1] == idx_to_del):
                    edges_to_del.add(edge)
            for edge in edges_to_del:
                self._edges.remove(edge)
            self._cleanup()

    def get_edges(self) -> set:
        """Returns the set of edges within this graph."""
        return self._edges

    def set_edges(self, new_edges: set[tuple[int]]) -> None:
        """Sets this graph's edge set to the given set."""
        self._edges = new_edges
        self._cleanup()

    def add_edge(self, new_edge: tuple[int]) -> None:
        """Adds a given edge to this graph."""
        self._edges.add(new_edge)
        self._cleanup()

    def del_edge(self, edge_to_del: tuple[int]) -> None:
        """Deletes a given edge within this graph."""
        if edge_to_del in self._edges:
            self._edges.remove(edge_to_del)
            self._cleanup()

    def check_edge(self, edge_to_check: tuple[int]) -> None:
        """Checks if a given edge is within this graph."""
        return edge_to_check in self._edges

    def get_rep(self) -> int:
        """Returns the index of this graph's representative vertex."""
        return self._rep

    def set_rep(self, idx: int) -> None:
        """Sets this graph's representative vertex to the object at the given
        index. Raises a ValueError if that object is not a child of this
        graph."""
        if idx in self.vertex_indices():
            self._rep = idx
        else:
            raise ValueError(
                f"The object at index {idx} is not a child of this graph."
            )  # noqa

    def get_adjacent(self, idx: int) -> set:
        """Returns a set containing the indices of all vertices adjacent to the
        vertex at the given index."""
        out = set([])
        if idx in self.vertex_indices():
            for edge in self._edges:
                if edge[0] == idx:
                    out.add(edge[1])
                elif edge[1] == idx:
                    out.add(edge[0])
        return out
