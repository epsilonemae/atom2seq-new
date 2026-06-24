from atom2seq.classes import Mol


def verify_mol(molecule: Mol) -> tuple[bool | str]:
    if len(molecule.get_atoms()) < 10:
        return (
            False,
            f"it only has {len(molecule.get_atoms())} atoms, and "
            "is too small to be a peptide.",
        )
    num_elements = {}
    req_atoms = ["H", "C", "O", "N"]
    valid_atoms = ["H", "C", "O", "N", "S"]
    for atom in molecule.get_atoms():
        if atom.symbol not in valid_atoms:
            return (
                False,
                f"{atom.symbol} is not an atom found in any encoded amino acid.",  # noqa
            )
        if atom.symbol in num_elements.keys():
            num_elements[atom.symbol] += 1
        else:
            num_elements[atom.symbol] = 1
    for elt in req_atoms:
        if elt not in num_elements.keys():
            return (
                False,
                f"it has no {elt} atoms in it, so it cannot be a " "peptide.",
            )
    return (True, "This is a possible peptide.")
