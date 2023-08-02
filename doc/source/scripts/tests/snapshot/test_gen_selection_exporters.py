from gen_selection_exporters import SelectionExporterWriter
from MDAnalysis.core import selection as sel
from unittest.mock import patch

def test_SelectionExporterWriter(snapshot):
    with patch("builtins.open"):
        se = SelectionExporterWriter()
    assert se.lines == snapshot
