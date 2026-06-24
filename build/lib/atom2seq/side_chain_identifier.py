from atom2seq.backbone_identifier import find_atom
from atom2seq.classes import Atom, Mol

aa_dict = {
    # (no. of C*, no. of H*, no. of N*, no. of O*, no. of S*): symbol
    # *in side chain
    (1, 3, 0, 0, 0): "A",
    (1, 3, 0, 0, 1): "C",
    (2, 3, 0, 2, 0): "D",
    (3, 5, 0, 2, 0): "E",
    (7, 7, 0, 0, 0): "F",
    (0, 1, 0, 0, 0): "G",
    (4, 6, 2, 0, 0): "H",
    # Leucine and Isoleucine are not distinguishable by number of atoms alone
    (4, 9, 0, 0, 0): "I/L",
    (4, 10, 1, 0, 0): "K",
    (3, 7, 0, 0, 1): "M",
    (2, 4, 1, 1, 0): "N",
    (3, 5, 0, 0, 0): "P",
    (3, 6, 1, 1, 0): "Q",
    (4, 12, 3, 0, 0): "R",
    (1, 3, 0, 1, 0): "S",
    (2, 5, 0, 1, 0): "T",
    (3, 7, 0, 0, 0): "V",
    (9, 8, 1, 0, 0): "W",
    (7, 7, 0, 1, 0): "Y",
}


def update_idx_del(molecule: Mol, current_idx: int, del_sym: str) -> int:
    to_del = []
    old_idx = False
    new_current_idx = False
    for i in molecule.get_bonded(current_idx):
        if molecule.get_backbone()[i]:
            old_idx = current_idx
            new_current_idx = i
        if molecule.get_atoms()[i].symbol == del_sym:
            to_del.append(i)
    for i in sorted(to_del, reverse=True):
        print(f"i in to_del: {molecule.get_atoms()[i]} is going to be deleted")
        molecule.del_atom(i)
        if old_idx > i:
            old_idx -= 1
        if new_current_idx > i:
            new_current_idx -= 1
    print(f"old_idx: {molecule.get_atoms()[old_idx]} is going to be deleted")
    molecule.del_atom(old_idx)
    if new_current_idx > old_idx:
        current_idx = new_current_idx - 1
    else:
        current_idx = new_current_idx
    return current_idx


def update_idx_delH_find_Rgroup(molecule: Mol, current_idx: int) -> tuple[int]:
    to_del = []
    old_idx = False
    new_current_idx = False
    r_group_idx = False
    for i in molecule.get_bonded(current_idx):
        if molecule.get_backbone()[i]:
            old_idx = current_idx
            new_current_idx = i
        # The only two things bonded to the alpha carbon are a carbon with 3
        # bonds (in the carboxylic acid) and the start of the R-group.
        elif len(molecule.get_bonded(i)) != 3:
            print("Found the R-group!")
            r_group_idx = i
        elif molecule.get_atoms()[i].symbol == "H":
            to_del.append(i)
    for i in sorted(to_del, reverse=True):
        print(f"i in to_del: {molecule.get_atoms()[i]} is going to be deleted")
        molecule.del_atom(i)
        if old_idx > i:
            old_idx -= 1
        if new_current_idx > i:
            new_current_idx -= 1
    print(f"old_idx: {molecule.get_atoms()[old_idx]} is going to be deleted")
    molecule.del_atom(old_idx)
    if new_current_idx > old_idx:
        current_idx = new_current_idx - 1
    else:
        current_idx = new_current_idx
    return (current_idx, r_group_idx)


def ile_or_leu(r_group):
    out = "I/L"
    for i in range(len(r_group)):
        atom = r_group.get_atoms()[i]
        # Checks how many carbons each carbon is bonded to. If one of
        # them has 3, then it is leucine.
        if atom.symbol == "C":
            bonded_carbons = find_atom(r_group, [i], "C")
            if len(bonded_carbons) == 3:
                out = "L"
    if out == "I/L":
        out = "I"
    return out


def id_side_chains(molecule: Mol) -> list[str]:
    print(f"{molecule=}")
    current_idx = molecule.get_n_term()
    out = []
    num_iter = 0
    # Makes a list containing only the atoms in the backbone.
    backbone_atoms = []
    for i in range(len(molecule.get_atoms())):
        if molecule.get_backbone()[i]:
            backbone_atoms.append(molecule.get_atoms()[i])
    # Makes a list of only the carbons in the backbone.
    backbone_carbons = [atom for atom in backbone_atoms if atom.symbol == "C"]
    # Alpha carbons make up exactly half of the carbons in the backbone.
    num_alpha_carbons = int(len(backbone_carbons) / 2)
    # The number of sidechains is the same as the number of alpha carbons.
    for _ in range(num_alpha_carbons):
        num_iter += 1
        print(f"{num_iter=}, {len(molecule.get_atoms())=}")
        # Finds the next step in the backbone (it will be the alpha carbon) and
        # deletes this nitrogen as well as any hydrogens connected to it.
        current_idx = update_idx_del(molecule, current_idx, "H")
        # Finds how many hydrogens are bonded to the alpha carbon. If it is 2,
        # this is a glycine and the R-group will be deleted, so we set the
        # R-group prematurely.
        bonded_syms = []
        for i in molecule.get_bonded(current_idx):
            bonded_syms.append(molecule.get_atoms()[i].symbol)
        numH = bonded_syms.count("H")
        r_group = False
        r_group_idx = False
        if numH == 2:
            print("It's a glycine!")
            current_idx = update_idx_del(molecule, current_idx, "H")
            r_group = Mol([Atom("H", (0, 0, 0))], [])
            print(f"{r_group=}")
        else:
            print("It's not a glycine!")
            # Deletes the alpha carbon and all hydrogens bonded to it, and
            # extracts the index of the r group.
            current_idx, r_group_idx = update_idx_delH_find_Rgroup(
                molecule, current_idx
            )
            print(f"{r_group_idx=}")
        print(f"{r_group=}, {r_group_idx=}")
        if (not r_group) and (r_group_idx):
            print(f"{r_group_idx=}, {molecule=}")
            r_group_idcs = molecule.find_submol(r_group_idx)
            r_group = molecule.split_submol(r_group_idx)
            for idx in r_group_idcs:
                if current_idx > idx:
                    current_idx -= 1
        r_group_makeup = [0, 0, 0, 0, 0]
        for atom in r_group.get_atoms():
            if atom.symbol == "C":
                r_group_makeup[0] += 1
            elif atom.symbol == "H":
                r_group_makeup[1] += 1
            elif atom.symbol == "N":
                r_group_makeup[2] += 1
            elif atom.symbol == "O":
                r_group_makeup[3] += 1
            elif atom.symbol == "S":
                r_group_makeup[4] += 1
        symbol = aa_dict[tuple(r_group_makeup)]
        if symbol == "I/L":
            symbol = ile_or_leu(r_group)
        out.append(symbol)
        current_idx = update_idx_del(molecule, current_idx, "O")
    return out
