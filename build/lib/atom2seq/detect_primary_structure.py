from atom2seq import parsers
from atom2seq.backbone_identifier import label_backbone
from atom2seq.bond_detector import bond_mol
from atom2seq.side_chain_identifier import id_side_chains
from atom2seq.verifier import verify_mol


def detect_primary_structure(filename: str, filetype: str):
    """Given a file and its type, returns the primary structure of the protein
    whose information is stored in that file."""
    # Reads the file using the parser for the filetype passed in. Allows a .xyz
    # filetype to be passed in as "xyz", "XYZ", ".xyz", or ".XYZ"
    molecule = False
    filetype.strip(".")
    if filetype.lower() == "gjf":
        molecule = parsers.parse_gjf(filename)
    elif filetype.lower() == "xyz":
        molecule = parsers.parse_xyz(filename)
    elif filetype.lower() == "pdb":
        molecule = parsers.parse_pdb(filename)
    elif filetype.lower() == "cif":
        molecule = parsers.parse_cif(filename)
    elif filetype.lower() == "nwc":
        molecule = parsers.parse_cif(filename)
    else:
        raise ValueError(f"{filetype.lower()} is not a supported filetype.")

    # Verifies that this is a potential peptide, raising a ValueError if it is
    # not.
    is_valid, text = verify_mol(molecule)
    if not is_valid:
        raise ValueError(
            f"The molecule in {filename} is not a valid protein, as {text}"
        )

    # Adds bonds, a backbone, and an N-terminus to the molecule.
    bond_mol(molecule)
    label_backbone(molecule)

    # Identifies the sidechains.
    side_chain_list = id_side_chains(molecule)

    # Formats the list of amino acids as a string.
    out = ""
    for symbol in side_chain_list:
        out += symbol

    return out
