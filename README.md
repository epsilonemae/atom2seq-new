# atom2seq

Detects primary structure of proteins given the Cartesian coordinates of the 
atoms.

# Installation

Installation instructions for users are coming soon!

## Development
All commands assume a Linux-like command line.

1. Obtain the repo from github
```
git clone https://github.com/SIMCODES-ISU/atom2seq.git
```

2. Create a Python virtual environment
```
cd atom2seq
python3 -m venv .venv
```

3. Activate the virtual environment
```
source .venv/bin/activate
```

4. Install the development dependencies
```
pip install ".[dev]"
```

5. Editable install the repo
```
pip install -e .
```

You can then test the repo by running
```
pytest
```

# Repository Table of Contents
> `docs`
> `src`
    > `atom2seq`
        > `__init__.py` added by Ryan M. Richard
        > `atom_class.py` added by Aspen A.M. Meissner
        > `connect_groups.py` added by Aspen A.M. Meissner
        > `connectivity_table_class.py` added by Aspen A.M. Meissner
        > `file2seq.py` added by Aspen A.M. Meissner
        > `get_pseq.py` added by Aspen A.M. Meissner
        > `group_class.py` added by Aspen A.M. Meissner
        > `grouper.py` added by Aspen A.M. Meissner
        > `mol_class.py` added by Aspen A.M. Meissner
        > `parsers.py` added by Aspen A.M. Meissner
        > `rgroup_class.py` added by Aspen A.M. Meissner
> `tests`
    > `assets`
        > `parser_tests`
            > `water_extra_info.cif` added by Aspen A.M. Meissner
            > `water_extra_info.pdb` added by Aspen A.M. Meissner
            > `water_extra_lines.xyz` added by Aspen A.M. Meissner
            > `water_no_lines.xyz` added by Aspen A.M. Meissner
            > `water_no_number_of_atoms.xyz` added by Aspen A.M. Meissner
            > `water.cif` added by Aspen A.M. Meissner
            > `water.gjf` added by Aspen A.M. Meissner
            > `water.pdb` added by Aspen A.M. Meissner
            > `water.xyz` added by Aspen A.M. Meissner
        > `glycine_optimized.xyz` added by Ella Bushman
        > `triY.nwc` added by Ella Bushman
        > `ubiquitin_full.py` added by Aspen A.M. Meissner
    > `test_atom_class.py` added by Aspen A.M. Meissner
    > `test_connect_groups.py` added by Aspen A.M. Meissner
    > `test_connectivity_table_class.py` added by Aspen A.M. Meissner
    > `test_file2seq.py` added by Aspen A.M. Meissner
    > `test_get_pseq.py` added by Aspen A.M. Meissner
    > `test_group_class.py` added by Aspen A.M. Meissner
    > `test_grouper.py` added by Aspen A.M. Meissner
    > `test_mol_class.py` added by Aspen A.M. Meissner
    > `test_parsers.py` added by Aspen A.M. Meissner
    > `test_rgroup_class.py` added by Aspen A.M. Meissner
    > `test_rgroup_symbol.py` added by Aspen A.M. Meissner
> `.gitignore` added by Ryan M. Richard
> `LICENSE` added by Ryan M. Richard
> `pyproject.toml` added by Ryan M. Richard
> `README.md` added by Ryan M. Richard

# To-Do
> Make `bond_mol` in `parsers.py` more efficient
> Make a PDB/CIF cleanup tool (currently these cannot be run due to extraneous information like floating water molecules)

# Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant No. 2348724.