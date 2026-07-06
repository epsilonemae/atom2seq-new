from atom2seq.fxnal_group_class import FxnalGroup
from atom2seq.indexed_graph_class import IndexedGraph


class Mol(IndexedGraph):
    def __init__(self, groups, bonds, parent: int = -1, idx: int = -1):
        super().__init__(groups, bonds, parent, idx)
        self._backbone = set()
        self._nterm = -1
