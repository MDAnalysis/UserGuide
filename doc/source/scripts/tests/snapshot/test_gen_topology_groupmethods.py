from unittest.mock import patch

from gen_topology_groupmethods import TransplantedMethods


def test_TransplantedMethods(snapshot):
    with patch("builtins.open"):
        tm = TransplantedMethods()
    assert tm.table_writer.lines == snapshot
