from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.rgroup_class import RGroup


# This is being done with a class to avoid having numerous functions that take
# in the backbone and aas lists as inputs.
class Grouped:
    def __init__(self, groups: set[Group], bonds: ConnectivityTable):
        self.groups = sorted(list(groups))
        for i in range(len(self.groups)):
            self.groups[i].set_idx(i)
        self.bonds = bonds
        self.backbone = []
        self.aas = []

    def group_symbols(self, idx: int) -> list:
        """Returns a sorted list containing the symbols of all atoms in the
        group at a given index."""
        return sorted([atom.symbol for atom in self.groups[idx].atom_list()])

    def id_nterms(self) -> list:
        """Returns a list of the indices of all NH2s and NHs within this
        object."""
        potential_nterms = []
        # Loops over all groups and extracts any NH2 or NH.
        for group in self.groups:
            if self.group_symbols(group.get_idx()) == ["H", "H", "N"]:
                potential_nterms.append(group.get_idx())
            elif self.group_symbols(group.get_idx()) == ["H", "N"]:
                potential_nterms.append(group.get_idx())
        return potential_nterms

    def bb_iter(self, current_idx: int) -> list | None:
        """Recursively finds each amino acid's side chain."""
        ch_step = False
        rgroup_idx = False
        self.backbone.append(current_idx)
        # if we are in proline, finds CH
        if self.group_symbols(current_idx) == ["N"]:
            for gp in self.bonds.get_paired(current_idx):
                if self.group_symbols(gp) == ["C", "H"]:
                    if gp not in self.backbone:
                        ch_step = gp
                        self.backbone.append(ch_step)
                    rgroup_idx = -1
        # finds CH or CH2 (glycine)
        else:
            for gp in self.bonds.get_paired(current_idx):
                if {"C", "H"}.issubset(set(self.group_symbols(gp))):
                    if gp not in self.backbone:
                        ch_step = gp
                        self.backbone.append(ch_step)
                    rgroup_idx = -1
        if ch_step:
            to_call = False
            to_return = False
            # look for next step, noting the r-group idx as well
            for gp in self.bonds.get_paired(ch_step):
                if (gp != current_idx) and (gp not in self.backbone):
                    # If this doesn't contain a carbonyl, it's the R-group
                    if not {"C", "O"}.issubset(set(self.group_symbols(gp))):
                        rgroup_idx = gp
                    else:
                        # if this is a cooh, we are done
                        if self.group_symbols(gp) == ["C", "H", "O", "O"]:
                            to_return = True
                            self.backbone.append(gp)
                        else:
                            if "N" not in self.group_symbols(gp):
                                # the only valid case where we are here is if
                                # the next aa is proline.
                                self.backbone.append(gp)
                                for new_gp in self.bonds.get_paired(gp):
                                    if self.group_symbols(new_gp) == ["N"]:
                                        to_call = new_gp
                            else:
                                # now the only case remaining is that the next
                                # group is Amide
                                to_call = gp
            if type(to_call) is int:
                self.rgroup_ider(rgroup_idx)
                return self.bb_iter(to_call)
            elif to_return:
                self.rgroup_ider(rgroup_idx)
                return self.aas

    def rgroup_ider(self, rgroup_idx):
        """Identifies an R-group based on the index of its beta carbon's
        group."""
        # Checks the Glycine edge case
        if rgroup_idx == -1:
            self.aas.append("G")
        else:
            # Initializes the RGroup object and adds the beta carbon's group to
            # it
            rgroup = RGroup(set(), ConnectivityTable(set()))
            rgroup.add_group(self.groups[rgroup_idx])
            complete = False
            # Loops until we have found the entire side chain
            while not complete:
                # Notes down the groups currently in the RGroup
                old_rgroup_groups = rgroup.group_list()
                # Loops over each group in the RGroup, then over each of the
                # groups bonded to that group. Checks if that group is in the
                # backbone - if it isn't, adds it to the RGroup and bonds it
                # within the RGroup to the group it is bonded to.
                to_add = []
                for group in rgroup.get_groups():
                    for gp in self.bonds.get_paired(group.get_idx()):
                        if gp not in self.backbone:
                            to_add.append(self.groups[gp])
                            rgroup.get_bonds().add_pair((group.get_idx(), gp))
                for group in to_add:
                    rgroup.add_group(group)
                # If the groups at the end of this loop are the same as the
                # groups at the beginning, we can exit the while loop.
                if rgroup.group_list() == old_rgroup_groups:
                    complete = True
            self.aas.append(rgroup.symbol())


def get_pseq(groups: set[Group], bonds: ConnectivityTable):
    molecule = Grouped(groups, bonds)
    nterms = molecule.id_nterms()
    if len(nterms) == 0:
        raise ValueError(
            "This is not a valid protein, as it has no potential N-termini."
        )
    found = False
    for nterm in nterms:
        # This will return None if it couldn't find a backbone off of this
        # N-term
        found = molecule.bb_iter(nterm)
        if found:
            return found
        molecule.backbone = []
        molecule.aas = []
    # This ValueError will only be raised if none of the N-terms found a
    # backbone.
    raise ValueError(
        "This is not a valid protein, as the backbone could not be found."
    )  # noqa
