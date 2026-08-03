from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group


def connect_groups(
    groups: set[Group], bonds: ConnectivityTable
) -> ConnectivityTable:  # noqa
    """Takes in a set of groups and the bonds between the atoms split between
    these groups to produce a ConnectivityTable of the bonds between the
    groups."""
    # Starts with an empty ConnectivityTable
    out = ConnectivityTable(set())
    # Assigns an index to each group
    groups = sorted(list(groups))
    for i in range(len(groups)):
        groups[i].set_idx(i)
    # Creates a dictionary from atom index to index of the group it is in
    in_what_group = {}
    for group in groups:
        for atom in group.atom_list():
            in_what_group[atom.get_idx()] = group.get_idx()
    # Loops over all bonds and adds a bond between the groups if they are not
    # the same
    for bond in bonds.get_pairs():
        if in_what_group[bond[0]] != in_what_group[bond[1]]:
            out.add_pair((in_what_group[bond[0]], in_what_group[bond[1]]))
    return out
