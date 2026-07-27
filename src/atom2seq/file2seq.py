from atom2seq.connect_groups import connect_groups
from atom2seq.get_pseq import get_pseq
from atom2seq.grouper import group_mol
from atom2seq.parsers import parse_cif, parse_gjf, parse_pdb, parse_xyz


def file2seq(filename: str, filetype: str):
    print("running file2seq")
    molecule = False
    # The filetype checker is being done this way so that it can be very easily
    # expanded in the future, either by SIMCODES or by the end user.
    supported_types = {
        "CIF": parse_cif,
        "GJF": parse_gjf,
        "PDB": parse_pdb,
        "XYZ": parse_xyz,
    }
    if filetype.upper() in supported_types:
        print("found supported type")
        molecule = supported_types[filetype.upper()](filename)
        print("parsed")
    else:
        raise ValueError(
            f"The .{filetype.upper()} filetype is not yet supported. Try one "
            f"of the following filetypes: {(supported_types.keys())}"
        )
    groups = group_mol(molecule)
    print("grouped")
    group_bonds = connect_groups(groups, molecule.get_bonds())
    print("groups bonded")
    return get_pseq(groups, group_bonds)
