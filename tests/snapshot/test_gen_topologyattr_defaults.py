
from doc.source.scripts.gen_topologyattr_defaults import TopologyDefaults
from MDAnalysis.core import selection as sel


def test_TopologyDefaults(snapshot):
    td = TopologyDefaults()
    assert td.lines == snapshot
