from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


class Backbone:
    """A class representing a functional group. Supports checking equality."""

    def __init__(self, groups: set[Group], bonds: ConnectivityTable):
        self._groups = groups
        self._bonds = bonds
        for i in range(len(self._groups)):
            self._groups[i].set_idx(i)

    def __repr__(self):
        return f"Mol({self._groups}, {self._bonds})"

    def __eq__(self, other):
        return (self._groups == other.get_groups()) and (
            self._bonds == other.get_bonds()
        )  # noqa

    def __hash__(self):
        return hash((self._groups, self._bonds))

    def get_bonds(self) -> ConnectivityTable:
        """Returns the ConnectivityTable of bonds."""
        return self._bonds

    def get_groups(self) -> list[Group]:
        """Returns the list of groups."""
        return self._groups
