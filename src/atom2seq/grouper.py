from atom2seq.mol_class import Mol


class GroupedMol(Mol):
    def __init__(self, molecule):
        super().__init__(molecule.get_atoms(), molecule.get_bonds())
        self.groups = set()
        self.grouped = set()

    def auto_group(self):
        for atom in self.atom_list():
            # Need to check individually each time so that only one ever gets
            # run on an atom. Note self.grouped is getting updated by the
            # detect functions.
            if atom.get_idx() not in self.grouped:
                self.detectAmd(atom.get_idx())
            if atom.get_idx() not in self.grouped:
                self.detectCOOH(atom.get_idx())
            if atom.get_idx() not in self.grouped:
                self.detectCO(atom.get_idx())
            if atom.get_idx() not in self.grouped:
                self.detectXHn(atom.get_idx())

    def detectAmd(self, idx: int, is_initial: bool = False):
        if is_initial:
            to_group = [idx]
            carbon = list(self._bonds.get_paired(idx))[0]  # only one thing
            to_group.append(carbon)
            for i in self._bonds.get_paired(carbon):
                if self.get_atom(i).symbol == "N":
                    to_group.append(i)
                    for j in self._bonds.get_paired(i):
                        if self.get_atom(j).symbol == "H":
                            to_group.append(j)
            if (len(to_group) == 4) or (len(to_group) == 5):
                for i in to_group:
                    self.grouped.add(i)
                self.groups.add(self.group_atoms(to_group))
        else:
            atom = self.get_atom(idx)
            if atom.symbol == "H":
                print(atom)
                print(self._bonds.get_paired(idx))
                print(self._bonds)
                bonded = list(self._bonds.get_paired(idx))[0]
                self.detectAmd(bonded)
            if atom.symbol == "O":
                if len(self._bonds.get_paired(idx)) == 1:
                    self.detectAmd(idx, True)
            elif atom.symbol == "C":
                for i in self._bonds.get_paired(idx):
                    if self.get_atom(i).symbol == "O":
                        self.detectAmd(i)
            elif atom.symbol == "N":
                for i in self._bonds.get_paired(idx):
                    if self.get_atom(i).symbol == "C":
                        if len(self._bonds.get_paired(i)) == 3:
                            self.detectAmd(i)

    def detectCOOH(self, idx: int, is_initial: bool = False):
        if is_initial:
            to_group = [idx]
            carbon = list(self._bonds.get_paired(idx))[0]  # only bonded to one
            to_group.append(carbon)
            for i in self._bonds.get_paired(carbon):
                if self.get_atom(i).symbol == "O":
                    if i not in to_group:
                        to_group.append(i)
                    for j in self._bonds.get_paired(i):
                        if self.get_atom(j).symbol == "H":
                            to_group.append(j)
            if len(to_group) == 4:
                for i in to_group:
                    self.grouped.add(i)
                self.groups.add(self.group_atoms(to_group))
        else:
            atom = self.get_atom(idx)
            if atom.symbol == "H":
                bonded = list(self._bonds.get_paired(idx))[0]
                self.detectCOOH(bonded)
            if atom.symbol == "O":
                if len(self._bonds.get_paired(idx)) == 1:
                    self.detectCOOH(idx, True)
                elif len(self._bonds.get_paired(idx)) == 2:
                    for i in self._bonds.get_paired(idx):
                        if self.get_atom(i).symbol == "C":
                            self.detectCOOH(i)
            elif atom.symbol == "C":
                for i in self._bonds.get_paired(idx):
                    if self.get_atom(i).symbol == "O":
                        if len(self._bonds.get_paired(i)) == 1:
                            self.detectCOOH(i, True)

    def detectCO(self, idx: int):
        atom = self.get_atom(idx)
        if atom.symbol == "O":
            if len(self._bonds.get_paired(idx)) == 1:
                i = list(self._bonds.get_paired(idx))[0]
                if self.get_atom(i).symbol == "C":
                    to_group = [idx, i]
                    for i in to_group:
                        self.grouped.add(i)
                    self.groups.add(self.group_atoms(to_group))
        elif atom.symbol == "C":
            for i in self._bonds.get_paired(idx):
                new_atom = self.get_atom(i)
                if new_atom.symbol == "O":
                    self.detectCO(i)

    def detectXHn(self, idx: int):
        atom = self.get_atom(idx)
        if atom.symbol == "H":
            center = list(self._bonds.get_paired(idx))[0]
            if center in self.grouped:
                self.groups.add(self.group_atoms([idx]))
            else:
                self.detectXHn(center)
        else:
            paired = list(self._bonds.get_paired(idx))
            to_group = []
            for i in paired:
                if self.get_atom(i).symbol == "H":
                    to_group.append(i)
            to_group.append(idx)
            for i in to_group:
                self.grouped.add(i)
            self.groups.add(self.group_atoms(to_group))


def group_mol(molecule):
    gmol = GroupedMol(molecule)
    gmol.auto_group()
    return gmol.groups
