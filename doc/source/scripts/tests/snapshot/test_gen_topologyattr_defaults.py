from unittest.mock import patch

from gen_topologyattr_defaults import TopologyDefaults
from MDAnalysis.core import selection as sel


def test_TopologyDefaults(snapshot):
    with patch("builtins.open"):
        td = TopologyDefaults()
    assert td.lines == snapshot
