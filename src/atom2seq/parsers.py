from atom2seq.atom_class import Atom
from atom2seq.connectivity_table_class import ConnectivityTable
from atom2seq.mol_class import Mol


def file_base(filename: str) -> list[str]:
    """Returns the stripped .readlines() of a passed file."""
    file = open(filename, "r")
    contents = file.readlines()
    file.close()
    contents = [line.strip() for line in contents]
    return contents


def parser_base(contents: list[list[str]]) -> Mol:
    """Takes in a list of lists of atomic symbols and their coordinates and
    returns a Mol containing all those atoms."""
    # Removes blank lines
    contents = [line for line in contents if line]
    new_contents = []
    for listy in contents:
        to_append = []
        for elt in listy:
            # If this elt in a float, changes it to be a float
            if elt.replace(".", "").replace("-", "").isdigit():
                to_append.append(float(elt))
            else:
                to_append.append(elt)
        new_contents.append(to_append)
    contents = new_contents
    # Creates a molecule containing each of these atoms and bonds it.
    atoms = set([Atom(listy[0], tuple(listy[1:])) for listy in contents])
    molecule = Mol(atoms, ConnectivityTable(set()))
    bond_mol(molecule)
    return molecule


def bond_mol(molecule: Mol) -> None:
    """Takes in an unbonded Mol and bonds it according to the covalent radii of
    atoms."""
    # Doing it by looping over all pairs instead of a KDTree because I simply
    # could not get the KDTree to work with any reasonably small k. It's
    # possible we could write our own with a variable search depth as opposed
    # to fixed.
    radii = {"H": 0.31, "O": 0.66, "N": 0.71, "C": 0.76, "S": 1.05}
    max_bonds = {"H": 1, "O": 2, "N": 3, "C": 4, "S": 2}
    # Loops over every pair of atoms
    for atom1 in molecule.atom_list():
        sym1 = atom1.symbol
        for atom2 in molecule.atom_list():
            sym2 = atom2.symbol
            # If the first atom is already bonded to the max number of atoms it
            # can be bonded to, skips it
            bonded1 = molecule.get_bonds().get_paired(atom1.get_idx())
            if len(bonded1) == max_bonds[atom1.symbol]:
                continue
            # If the distance between them is less than 1.15 times the sum of
            # their covalent radii, bonds them
            max_dist = 1.15 * (radii[sym1] + radii[sym2])
            if (
                0.5
                <= molecule.dist(atom1.get_idx(), atom2.get_idx())
                <= max_dist  # noqa
            ):  # noqa
                molecule.get_bonds().add_pair(
                    (atom1.get_idx(), atom2.get_idx())
                )  # noqa


def parse_gjf(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .gjf format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    # GJFs are not split by line so we have to reconstruct it
    new_contents = ""
    for line in contents:
        new_contents += line
    contents = new_contents
    # Splits by \
    contents = contents.split("\\")
    # Splits each line by ,
    contents = [line.split(",") for line in contents]

    return parser_base(contents)


def parse_xyz(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .xyz format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    # Checking if the first line is the number of atoms. If it is, remove
    # the first line and any blank lines that come after it.
    if contents[0][0].isdigit():
        contents.pop(0)
    contents = [line.split() for line in contents]

    return parser_base(contents)


def parse_pdb(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .pdb format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    # Keeps each line if and only if it starts with ATOM
    contents = [line for line in contents if line[0:4] == "ATOM"]
    new_contents = []
    for line in contents:
        # Only keeps the A branch
        if line.split()[-8] == "A":
            # Extracts the symbol and coordinates
            new_contents.append([line.split()[-1], *line.split()[-6:-3]])

    return parser_base(new_contents)


def parse_cif(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .cif format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    contents = [line for line in contents if line[0:4] == "ATOM"]
    new_contents = []
    for line in contents:
        # Only keeps the A branch
        if line.split()[6] == "A":
            # Extracts the symbol and coordinates
            new_contents.append([line.split()[2], *line.split()[10:13]])

    return parser_base(new_contents)


def parse_nwc(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .nwc format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    contents = [line.split() for line in contents]
    new_contents = []
    for line in contents:
        # Extracts the symbol and coordinates
        new_contents.append([line[1], *line[3:]])

    return parser_base(new_contents)
