
from doc.source.scripts.gen_topologyparser_attrs import TopologyParsers, TopologyAttrs, ConnectivityAttrs
from MDAnalysis.core import selection as sel


def test_TopologyParsers(snapshot):
    top = TopologyParsers()
    assert top.lines == snapshot

def test_TopologyAttrs(snapshot):
    top = TopologyParsers()
    ta = TopologyAttrs(top.attrs)
    assert ta.lines == snapshot

def test_ConnectivityAttrs(snapshot):
    top = TopologyParsers()
    ca = ConnectivityAttrs(top.attrs)
    assert ca.lines == snapshot
