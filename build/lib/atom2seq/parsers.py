from atom2seq.classes import Atom, Mol


def file_base(filename: str) -> list[str]:
    """Returns the stripped .readlines() of a passed file."""
    file = open(filename, "r")
    contents = file.readlines()
    file.close()
    contents = [line.strip() for line in contents]
    return contents


def parser_base(contents: list[list[str | int]]) -> Mol:
    """Takes in a list of lists of atomic symbols and their coordinates and
    returns a Mol containing all those atoms."""
    # Removes blank lines, changes any integers to be ints, and then returns a
    # Mol containing atoms made from each element of contents and no bonds.
    contents = [line for line in contents if line]
    contents = [
        [
            (
                float(elt)
                if elt.replace(".", "").replace("-", "").isdigit()
                else elt  # noqa
            )  # noqa
            for elt in listy
        ]
        for listy in contents  # noqa
    ]
    atoms = [Atom(listy[0], tuple(listy[1:])) for listy in contents]
    return Mol(atoms, [])


def parse_gjf(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .gjf format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    new_contents = ""
    for line in contents:
        new_contents += line
    contents = new_contents
    contents = contents.split("\\")
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

    contents = [line for line in contents if line[0:4] == "ATOM"]
    contents = [
        [line.split()[-1], *line.split()[-6:-3]]
        for line in contents
        if line.split()[-8] == "A"
    ]

    return parser_base(contents)


def parse_cif(filename: str) -> Mol:
    """Parses the coordinates of a Molecule stored in .cif format and returns a
    Mol object containing those atoms."""
    contents = file_base(filename)

    contents = [line for line in contents if line[0:4] == "ATOM"]
    contents = [
        [line.split()[2], *line.split()[10:13]]
        for line in contents
        if line.split()[6] == "A"
    ]

    return parser_base(contents)


def parse_nwc(filename: str) -> Mol:
    contents = file_base(filename)

    contents = [[line.split()[1], *line.split()[3:]] for line in contents]

    return parser_base(contents)
