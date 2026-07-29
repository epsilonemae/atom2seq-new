from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.group_class import Group
from atom2seq.rgroup_class import RGroup


class Grouped:
    def __init__(self, groups: set[Group], bonds: ConnectivityTable):
        self.groups = sorted(list(groups))
        for i in range(len(self.groups)):
            self.groups[i].set_idx(i)
        self.bonds = bonds
        self.backbone = []
        self.aas = []

    def group_symbols(self, idx):
        return sorted([atom.symbol for atom in self.groups[idx].atom_list()])

    def id_nterms(self):
        potential_nterms = []
        for group in self.groups:
            if self.group_symbols(group.get_idx()) == ["H", "H", "N"]:
                potential_nterms.append(group.get_idx())
            elif self.group_symbols(group.get_idx()) == ["H", "N"]:
                potential_nterms.append(group.get_idx())
        return potential_nterms

    def bb_iter(self, current_idx: int):
        print(f"beginning of bb_iter: {self.bonds=}")
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
            if (to_call == 0) and (type(to_call) is int):
                print(self.backbone)
                self.rgroup_ider(rgroup_idx)
                return self.bb_iter(to_call)
            elif to_call:
                print(self.backbone)
                self.rgroup_ider(rgroup_idx)
                return self.bb_iter(to_call)
            elif to_return:
                print(self.backbone)
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
                print("==========================")
                print(rgroup)
                old_rgroup_groups = rgroup.group_list()
                to_add = []
                for group in rgroup.get_groups():
                    print(f"{group.get_idx()=}")
                    for gp in self.bonds.get_paired(group.get_idx()):
                        if gp not in self.backbone:
                            print(f"{gp=}")
                            to_add.append(self.groups[gp])
                            rgroup.get_bonds().add_pair((group.get_idx(), gp))
                for group in to_add:
                    rgroup.add_group(group)
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
    # This ValueError will only be raised if none of the N-terms found a
    # backbone.
    raise ValueError(
        "This is not a valid protein, as the backbone could not be found."
    )  # noqa
