from unittest.mock import patch

from gen_selection_exporters import SelectionExporterWriter


def test_SelectionExporterWriter(snapshot):
    with patch("builtins.open"):
        se = SelectionExporterWriter()
    assert se.table_writer.lines == snapshot
