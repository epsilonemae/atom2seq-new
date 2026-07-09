import warnings

from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


class Backbone:
    """A class representing a protein's backbone. Supports checking
    equality."""

    def __init__(self, groups: set[Group], bonds: ConnectivityTable):
        self._groups = groups
        self._bonds = bonds
        for i in range(len(self._groups)):
            self.group_list()[i].set_idx(i)

    def __repr__(self):
        return f"Backbone({self._groups}, {self._bonds})"

    def __eq__(self, other):
        return (self._groups == other.get_groups()) and (
            self._bonds == other.get_bonds()
        )  # noqa

    def group_list(self):
        return sorted(list(self._groups))

    def get_bonds(self) -> ConnectivityTable:
        """Returns the ConnectivityTable of bonds."""
        return self._bonds

    def get_groups(self) -> list[Group]:
        """Returns the list of groups."""
        return self._groups

    def add_group(self, new_group: Group) -> None:
        """Adds the given group to the backbone."""
        existing_idxs = {group.get_idx() for group in self._groups}
        new_idx = max(existing_idxs) + 1
        if new_group.get_idx() == -1:
            new_group.set_idx(new_idx)
        elif new_group.get_idx() in existing_idxs:
            warnings.warn(
                "The existing index is already in this backbone, so it is "
                f"being changed to {new_idx}."
            )
            new_group.set_idx(new_idx)
        self._groups.add(new_group)
