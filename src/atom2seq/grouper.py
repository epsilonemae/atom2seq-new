from atom2seq.mol_class import Mol


# Does this as a class so that we don't have functions that take in some list
# of atoms that have already been grouped.
class GroupedMol(Mol):
    def __init__(self, molecule):
        super().__init__(molecule.get_atoms(), molecule.get_bonds())
        self.groups = set()
        self.grouped = set()

    def auto_group(self):
        """Groups every atom into their respective functional groups."""
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
        """Detects an amide containing the given index. Adds the amide to the
        object's groups, adding nothing if the passed index's atom is not in an
        amide."""
        # If this is the atom that we identify the group from, identifies it
        # from there.
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
        # If this isn't, identifies a step along the path to that atom and
        # recursively calls the function.
        else:
            atom = self.get_atom(idx)
            if atom.symbol == "H":
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
        """Detects a COOH containing the given index. Adds the COOH to the
        object's groups, adding nothing if the passed index's atom is not in a
        COOH."""
        # If this is the atom that we identify the group from, identifies it
        # from there.
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
        # If this isn't, identifies a step along the path to that atom and
        # recursively calls the function.
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
        """Detects a carbonyl containing the given index. Adds the carbonyl to
        the object's groups, adding nothing if the passed index's atom is not
        in a carbonyl."""
        atom = self.get_atom(idx)
        # If this is the oxygen, identifies the whole group and moves on.
        if atom.symbol == "O":
            if len(self._bonds.get_paired(idx)) == 1:
                i = list(self._bonds.get_paired(idx))[0]
                if self.get_atom(i).symbol == "C":
                    to_group = [idx, i]
                    for i in to_group:
                        self.grouped.add(i)
                    self.groups.add(self.group_atoms(to_group))
        # If this is the carbon, recursively runs this function on the oxygen.
        elif atom.symbol == "C":
            for i in self._bonds.get_paired(idx):
                new_atom = self.get_atom(i)
                if new_atom.symbol == "O":
                    self.detectCO(i)

    def detectXHn(self, idx: int):
        """Detects a cluster containing the given index. Adds the cluster to
        the object's groups."""
        atom = self.get_atom(idx)
        # If this is a hydrogen, identifies the center.
        if atom.symbol == "H":
            center = list(self._bonds.get_paired(idx))[0]
            if center in self.grouped:
                self.groups.add(self.group_atoms([idx]))
            else:
                self.detectXHn(center)
        # Otherwise, finds all of the hydrogens bonded to this atom.
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
