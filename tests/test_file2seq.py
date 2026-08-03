from atom2seq.file2seq import file2seq

prefix = __file__.removesuffix("test_file2seq.py")


def test_glycine():
    gly = prefix + "assets/glycine_optimized.xyz"
    assert file2seq(gly, "XYZ") == "G"


def test_triY():
    triY = prefix + "assets/triY.nwc"
    assert file2seq(triY, "NWC") == "YYY"
