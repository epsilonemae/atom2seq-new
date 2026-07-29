from atom2seq.connect_groups import connect_groups
from atom2seq.file2seq import file2seq
from atom2seq.get_pseq import get_pseq
from atom2seq.group_class import Group
from atom2seq.grouper import group_mol
from atom2seq.ubiquitin import ubiquitin

print(
    file2seq(
        "/home/aspenamm/Documents/atom2seq-new/tests/assets/ubiquitin_full.xyz",  # noqa
        "XYZ",  # noqa
    )
)

# def symbol(group: Group):
#     symbols_dict = {}
#     for atom in group.atom_list():
#         if atom.symbol in symbols_dict:
#             symbols_dict[atom.symbol] += 1
#         else:
#             symbols_dict[atom.symbol] = 1
#     out = ""
#     for symbol in symbols_dict:
#         if symbols_dict[symbol] == 1:
#             out += symbol
#         else:
#             out += symbol + str(symbols_dict[symbol])
#     return out


# groups = group_mol(ubiquitin)
# bonds = connect_groups(groups, ubiquitin.get_bonds())
# pseq = get_pseq(groups, bonds)

# for letter in pseq:
#     print(letter, end="")
