
from doc.source.scripts.gen_unit_tables import write_unit_table
from unittest.mock import patch, mock_open
from pathlib import Path

def test_write_unit_table(snapshot):
    with patch.object(Path, "open", mock_open()) as mock_open_file:
        write_unit_table(filename=Path())

    lines = [args[0] for args in mock_open_file.return_value.write.call_args_list]
    assert lines == snapshot
