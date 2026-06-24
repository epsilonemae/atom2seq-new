from atom2seq.classes import Mol


def id_n_term(molecule: Mol) -> list[int]:
    out = []
    # loops over the atoms in the molecule.
    for i in range(len(molecule.get_atoms())):
        atom = molecule.get_atoms()[i]
        valid = False
        # only considers the atom if it is a nitrogen.
        if atom.symbol == "N":
            # makes list of the symbols of atoms bonded to nitrogen.
            bonds = molecule.get_bonded(i)
            bonded_syms = [molecule.get_atoms()[bond].symbol for bond in bonds]
            # initializes counters and helper variables.
            found_invalid = False
            h_counter = 0
            c_counter = 0
            for sym in bonded_syms:
                # Each hydrogen adds one to h_counter, carbon to c_counter, and
                # any other atom triggers the found_invalid flag.
                if sym == "H":
                    h_counter += 1
                elif sym == "C":
                    c_counter += 1
                else:
                    found_invalid = True
            # if the found_invalid flag was not triggered, then we check if it
            # is in a non-Proline AA.
            if not found_invalid:
                if c_counter == 1:
                    if h_counter == 2 or h_counter == 3:
                        valid = True
            if valid:
                out.append(i)
    print(f"N-terms = {out}")
    return out


def id_pro_n_term(molecule: Mol) -> list[int]:
    out = []
    # loops over the atoms in the molecule.
    for i in range(len(molecule.get_atoms())):
        atom = molecule.get_atoms()[i]
        valid = False
        # only considers the atom if it is a nitrogen.
        if atom.symbol == "N":
            # makes list of the symbols of atoms bonded to nitrogen.
            bonds = molecule.get_bonded(i)
            bonded_syms = [molecule.get_atoms()[bond].symbol for bond in bonds]
            # initializes counters and helper variables.
            found_invalid = False
            h_counter = 0
            c_counter = 0
            for sym in bonded_syms:
                # Each hydrogen adds one to h_counter, carbon to c_counter, and
                # any other atom triggers the found_invalid flag.
                if sym == "H":
                    h_counter += 1
                elif sym == "C":
                    c_counter += 1
                else:
                    found_invalid = True
            # if the found_invalid flag was not triggered, then we check if it
            # is in a Proline.
            if not found_invalid:
                j = 0
                if (c_counter == 2) and (h_counter == 1):
                    # finds all carbons bonded to the nitrogen
                    carbons = [
                        index
                        for index in molecule.get_bonded(i)
                        if molecule.get_atoms()[index].symbol() == "C"
                    ]
                    # finds a carbon that is bonded to the first carbon bonded
                    # to the nitrogen that is not bonded to an oxygen. there is
                    # only one of these in a real proline, so we can safely
                    # take the first one.
                    new_carbon = [
                        index
                        for index in molecule.get_bonded(carbons[0])
                        if (molecule.get_atoms()[index].symbol == "C")
                        and ("O" not in molecule.get_bonds(index))
                    ][0]
                    for idx1 in molecule.get_bonded(new_carbon):
                        for idx2 in molecule.get_bonded(carbons[1]):
                            if idx1 == idx2:
                                j += 1
                if j != 0:
                    valid = True
            if valid:
                out.append(i)
    print(f"Proline N-terms = {out}")
    return out


def find_atom(molecule: Mol, listy: list, sym: str) -> list[int]:
    idx = listy[-1]
    bonds = molecule.get_bonded(idx)
    out = []
    for i in bonds:
        if molecule.get_atoms()[i].symbol == sym:
            out.append(i)
    return out


def small_backbone_iter(
    molecule: Mol, listy: list, step: str
) -> list[list[int]]:  # noqa
    print(f"Finding a {step}")
    new_list = []
    for sublist in listy:
        next_atoms = find_atom(molecule, sublist, step)
        print(f"{sublist=}, {next_atoms=}")
        for atom in next_atoms:
            if atom not in sublist:
                if ((step == "O") and (len(next_atoms) > 1)) or (step != "O"):
                    print(f"index {atom} is a valid next step")
                    to_add = sublist + [atom]
                    new_list.append(to_add)
                else:
                    print(f"index {atom} is not a valid next step")
            else:
                print(f"index {atom} already in sublist")
    return new_list


def large_backbone_iter(
    molecule: Mol, listy: list
) -> tuple[list[list[int]], bool]:  # noqa
    plus_c = small_backbone_iter(molecule, listy, "C")
    if len(plus_c) != 0:
        plus_c = small_backbone_iter(molecule, plus_c, "C")
        if len(plus_c) != 0:
            plus_n = small_backbone_iter(molecule, plus_c, "N")
            if len(plus_n) != 0:
                return plus_n, False
            plus_o = small_backbone_iter(molecule, plus_c, "O")
            if len(plus_o) != 0:
                return plus_o, True
            else:
                raise ValueError(
                    "This molecule is not a protein as it has no backbone."
                )  # noqa
        else:
            raise ValueError(
                "This molecule is not a protein as it has no backbone."
            )  # noqa
    else:
        raise ValueError(
            "This molecule is not a protein as it has no backbone."
        )  # noqa


def find_backbone(molecule: Mol) -> list[int]:
    n_termini = id_n_term(molecule)
    paths = [[n_terminus] for n_terminus in n_termini]
    found = False
    i = 0
    while not found:
        paths, found = large_backbone_iter(molecule, paths)
        i += 1
        if i > 2 * len(molecule.get_atoms()):
            break
    if len(paths) == 0:
        n_termini = id_pro_n_term(molecule)
        paths = [[n_terminus] for n_terminus in n_termini]
        found = False
        i = 0
        while not found:
            paths, found = large_backbone_iter(molecule, paths)
            i += 1
            if i > 2 * len(molecule.get_atoms()):
                raise RecursionError("While loop maximum depth exceeded.")
        if len(paths) == 0:
            raise ValueError(
                "This molecule is not a protein as it has no backbone."
            )  # noqa
    return paths[0]


def label_backbone(molecule: Mol) -> None:
    backbone_idcs = find_backbone(molecule)
    molecule.set_backbone(backbone_idcs)
    molecule.set_n_term(backbone_idcs[0])
