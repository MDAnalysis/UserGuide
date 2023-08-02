
from doc.source.scripts.gen_unit_tables import write_unit_table
from unittest.mock import patch, mock_open

def test_write_unit_table(snapshot):
    with patch("builtins.open", new_callable=mock_open()) as mock_open_file:
        tables = write_unit_table(filename="fake.txt")

        raise ValueError(mock_open_file.return_value.call_args_list)
    assert tables == snapshot
