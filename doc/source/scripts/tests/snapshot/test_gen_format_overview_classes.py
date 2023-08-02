from gen_format_overview_classes import FormatOverview, FILE_TYPES, CoordinateReaders, SphinxClasses
from unittest.mock import patch

def test_FILE_TYPES():
    assert  FILE_TYPES.keys() >= {"DCD", "GRO", "XTC", "ITP"}


def test_FormatOverview(snapshot):
    with patch("builtins.open"):
        ov = FormatOverview()
    assert ov.lines == snapshot


def test_CoordinateReaders(snapshot):
    with patch("builtins.open"):
        cr = CoordinateReaders()
    assert cr.lines == snapshot


def test_SphinxClasses(snapshot):
    with patch("builtins.open"):
        sc = SphinxClasses("PDB")
    assert sc.lines == snapshot
