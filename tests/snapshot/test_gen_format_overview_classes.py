from doc.source.scripts.gen_format_overview_classes import FormatOverview, FILE_TYPES, CoordinateReaders, SphinxClasses


def test_FILE_TYPES():
    assert  FILE_TYPES.keys() >= {"DCD", "GRO", "XTC", "ITP"}

def test_FormatOverview(snapshot):
    ov = FormatOverview()
    assert ov.lines == snapshot


def test_CoordinateReaders(snapshot):
    cr = CoordinateReaders()
    assert cr.lines == snapshot


def test_SphinxClasses(snapshot):
    sc = SphinxClasses("PDB")
    assert sc.lines == snapshot
