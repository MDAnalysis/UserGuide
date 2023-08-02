from doc.source.scripts.gen_selection_exporters import SelectionExporterWriter
from MDAnalysis.core import selection as sel

def test_SelectionExporterWriter(snapshot):
    se = SelectionExporterWriter()
    assert se.lines == snapshot
