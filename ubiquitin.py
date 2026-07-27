from atom2seq.file2seq import file2seq

prefix = __file__.removesuffix("ubiquitin.py")

print(file2seq(prefix + "tests/assets/ubiquitin_full.xyz", "xyz"))
