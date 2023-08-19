from unittest.mock import patch

from gen_topologyattr_defaults import TopologyDefaults


def test_TopologyDefaults(snapshot):
    with patch("builtins.open"):
        td = TopologyDefaults()
    assert td.table_writer.lines == snapshot
