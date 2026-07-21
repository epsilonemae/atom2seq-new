from atom2seq.backbone_class import Backbone
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.rgroup_class import RGroup


class Grouped:
    def __init__(self, groups: set[Group], bonds: ConnectivityTable):
        self.groups = sorted(list(groups))
        self.bonds = bonds
        self.backbone = Backbone(set(), ConnectivityTable(set()))
        self.aas = []

    def group_symbols(self, idx):
        group = False
        for gp in self.groups:
            if gp.get_idx() == idx:
                group = gp
        return sorted([atom.symbol for atom in group.atom_list()])

    def id_nterms(self):
        potential_nterms = []
        for group in self.groups:
            if self.group_symbols(group.get_idx()) == ["H", "H", "N"]:
                potential_nterms.append(group)
            elif self.group_symbols(group.get_idx()) == ["H", "N"]:
                potential_nterms.append(group)
        return potential_nterms

    def bb_iter(self, current_idx: int):
        # finds CH or CH2 (glycine)
        rgroup_idx = False
        ch_step = False
        for gp in self.bonds.get_paired(current_idx):
            if {"C", "H"}.issubset(set(self.group_symbols(gp))):
                ch_step = gp
                rgroup_idx = -1
        if ch_step:
            to_call = False
            to_return = False
            # look for next step, noting the r-group idx as well
            for gp in self.bonds.get_paired(ch_step):
                if gp != current_idx:
                    # If this doesn't contain a carbonyl, it's the R-group
                    if not {"C", "O"}.issubset(set(self.group_symbols(gp))):
                        rgroup_idx = gp
                        # call some function that ids the rgroup
                    else:
                        if self.group_symbols(gp) == ["C", "H", "O", "O"]:
                            to_return = True
                        else:
                            # the only valid case where this is not true is if
                            # the next aa is proline.
                            if "N" not in self.group_symbols(gp):
                                for new_gp in self.bonds.get_paired(gp):
                                    if self.group_symbols(new_gp) == ["N"]:
                                        to_call = new_gp
                            else:
                                # now the only case remaining is Amide
                                to_call = new_gp
            if to_call:
                self.bb_iter(to_call)
                self.rgroup_ider(rgroup_idx)
            elif to_return:
                self.rgroup_ider(rgroup_idx)
                return self.aas

    def rgroup_ider(self, rgroup_idx):
        if rgroup_idx == -1:
            self.aas.append("G")
        else:
            rgroup = RGroup(set(), ConnectivityTable(set()))
            rgroup.add_group(self.groups[rgroup_idx])
            complete = False
            while not complete:
                for group in rgroup.get_groups():
                    for gp in self.bonds.get_paired(group.get_idx()):
                        check_gp = self.groups[gp]
                        if check_gp not in self.backbone:
                            rgroup.add_group(check_gp)


def id_nterms(groups: set[Group], bonds: ConnectivityTable):
    # The indices should already be set from the connect_groups() f'n.
    groups = sorted(list(groups))
    # Look for N-term (NH2 or NH in a Proline)
    potential_nterms = []
    for group in groups:
        symbols = [atom.symbol for atom in group.atom_list()]
        symbols = sorted(symbols)
        if (symbols == ["H", "H", "N"]) or (symbols == ["H", "N"]):
            potential_nterms.append(group)
    for nterm in potential_nterms:
        pass
        # loop over calling a function that:
        # finds CH
        # finds Amd
        # go from Amd
        # if no Amd: (we are either at C-term or in P)
        #   find C=O
        # if no C=O: (we are either at C-term or invalid)
        #   find COOH
        #   if no COOH:
        #       invalid
        #   find N
        #   if no N:
        #       invalid
        #   go from N
