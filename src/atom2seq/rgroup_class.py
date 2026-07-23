import warnings

from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


class RGroup:
    """A class representing an  R-group. Supports checking equality."""

    # The symbol_dict does not contain G since G is covered by a special case
    # in the backbone/AA ID'er.
    symbol_dict = {
        (("C", "H", "H", "H"),): "A",
        (("C", "H", "H"), ("H", "S")): "C",
        (("C", "H", "H"), ("C", "H", "O", "O")): "D",
        (("C", "H", "H"), ("C", "H", "H"), ("C", "H", "O", "O")): "E",
        (
            ("C",),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H", "H"),
        ): "F",
        (
            ("C",),
            ("C", "H"),
            ("C", "H"),
            ("C", "H", "H"),
            ("H", "N"),
            ("N",),
        ): "H",  # noqa
        (
            ("C", "H"),
            ("C", "H", "H"),
            ("C", "H", "H", "H"),
            ("C", "H", "H", "H"),
        ): "I/L",  # noqa
        (
            ("C", "H", "H"),
            ("C", "H", "H"),
            ("C", "H", "H"),
            ("C", "H", "H"),
            ("H", "H", "N"),
        ): "K",
        (("C", "H", "H"), ("C", "H", "H"), ("C", "H", "H", "H"), ("S",)): "M",
        (("C", "H", "H"), ("C", "H", "N", "O"), ("H",)): "N",
        (("C", "H", "H"), ("C", "H", "H"), ("C", "H", "H")): "P",
        (("C", "H", "H"), ("C", "H", "H"), ("C", "H", "N", "O"), ("H",)): "Q",
        (
            ("C",),
            ("C", "H", "H"),
            ("C", "H", "H"),
            ("C", "H", "H"),
            ("H", "H", "N"),
            ("H", "N"),
            ("H", "N"),
        ): "R",
        (("C", "H", "H"), ("H", "O")): "S",
        (("C", "H"), ("C", "H", "H", "H"), ("H", "O")): "T",
        (("C", "H"), ("C", "H", "H", "H"), ("C", "H", "H", "H")): "V",
        (
            ("C",),
            ("C",),
            ("C",),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H", "H"),
            ("H", "N"),
        ): "W",
        (
            ("C",),
            ("C",),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H"),
            ("C", "H", "H"),
            ("H", "O"),
        ): "Y",
    }

    def __init__(self, groups: set[Group], bonds: ConnectivityTable):
        self._groups = groups
        self._bonds = bonds
        for i in range(len(self._groups)):
            self.group_list()[i].set_idx(i)

    def __repr__(self):
        return f"RGroup({self.group_list()}, {self._bonds})"

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
        """Returns the set of groups."""
        return self._groups

    def symbol(self):
        key = []
        for group in self._groups:
            add_to_key = []
            for atom in group.get_atoms():
                add_to_key.append(atom.symbol)
            key.append(tuple(sorted(add_to_key)))
        symbol = self.symbol_dict[tuple(sorted(key))]
        if symbol == "I/L":
            for group in self._groups:
                group_symbols = [a.symbol for a in group.get_atoms()]
                if sorted(group_symbols) == ["C", "H"]:
                    idx = group.get_idx()
                    if len(self._bonds.get_paired(idx)) == 3:
                        return "L"
                    elif len(self._bonds.get_paired(idx)) == 2:
                        return "I"
                    else:
                        raise KeyError(
                            "This R-group is I/L, but is neither I nor L."
                        )  # noqa
        else:
            return symbol

    def add_group(self, new_group: Group) -> None:
        """Adds the given group to the R-Group."""
        existing_idxs = {group.get_idx() for group in self._groups}
        if existing_idxs:
            new_idx = max(existing_idxs) + 1
        else:
            new_idx = 0
        if new_group.get_idx() == -1:
            new_group.set_idx(new_idx)
        elif new_group.get_idx() in existing_idxs:
            new_group.set_idx(new_idx)
        self._groups.add(new_group)
