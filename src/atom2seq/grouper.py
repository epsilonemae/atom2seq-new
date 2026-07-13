from atom2seq.atom_class import Atom
from atom2seq.group_class import Group
from atom2seq.mol_class import Mol


class GroupedMol(Mol):
    def __init__(self, molecule):
        super().__init__(molecule.get_atoms(), molecule.get_bonds())
        self.groups = set()
        self.grouped = []

    def auto_group(self):
        for atom in self.atom_list():
            if atom.get_idx() not in self.grouped:
                self.detectInd(atom.get_idx())
                self.detectImd(atom.get_idx())
                self.detectPhOH(atom.get_idx())
                self.detectPh(atom.get_idx())
                self.detectAmd(atom.get_idx())
                self.detectCOOH(atom.get_idx())
                self.detectCO(atom.get_idx())
                self.detectXHn(atom.get_idx())

    def detectInd(self, idx: int):
        pass

    def detectImd(self, idx: int):
        pass

    def detectPhOH(self, idx: int):
        pass

    def detectPh(self, idx: int):
        pass

    def detectAmd(self, idx: int):
        pass

    def detectCOOH(self, idx: int):
        pass

    def detectCO(self, idx: int):
        pass

    def detectXHn(self, idx: int):
        if self.get_atom(idx).symbol == "H":
            center = self.get_bonds().get_paired(idx)[0]
            self.detectXHn(center)
        else:
            to_group = self.get_bonds().get_paired(idx).append(idx)
            for i in to_group:
                self.grouped.append(i)
            self.out.add(self.group_atoms(to_group))


def group_mol(molecule: Mol) -> set[Group]:
    out = set()
    grouped = []
    for atom in molecule.atom_list():
        if atom.get_idx() not in grouped:
            detectInd(molecule, atom.get_idx(), out)
            detectImd(molecule, atom.get_idx(), out)
            detectPhOH(molecule, atom.get_idx(), out)
            detectPh(molecule, atom.get_idx(), out)
            detectAmd(molecule, atom.get_idx(), out)
            detectCOOH(molecule, atom.get_idx(), out)
            detectCO(molecule, atom.get_idx(), out)
            detectXHn(molecule, atom.get_idx(), out)


def detectInd(molecule: Mol, idx: int, out: set, grouped: list):
    pass


def detectImd(molecule: Mol, idx: int, out: set, grouped: list):
    pass


def detectPhOH(molecule: Mol, idx: int, out: set, grouped: list):
    pass


def detectPh(molecule: Mol, idx: int, out: set, grouped: list):
    pass


def detectAmd(molecule: Mol, idx: int, out: set, grouped: list):
    pass


def detectCOOH(molecule: Mol, idx: int, out: set, grouped: list):
    pass


def detectCO(molecule: Mol, idx: int, out: set, grouped: list):
    if molecule.get_atom(idx).symbol == "C":
        for i in molecule.get_bonds().get_paired(idx):
            if molecule.get_atom(i).symbol == "O":
                detectCO(molecule, i, out)


def detectXHn(molecule: Mol, idx: int, out: set, grouped: list):
    if molecule.get_atom(idx).symbol == "H":
        center = molecule.get_bonds().get_paired(idx)[0]
        detectXHn(molecule, center, out)
    else:
        to_group = molecule.get_bonds().get_paired(idx).append(idx)
        out.add(molecule.group_atoms(to_group))
