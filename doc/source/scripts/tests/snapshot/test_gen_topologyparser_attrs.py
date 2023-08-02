
from gen_topologyparser_attrs import TopologyParsers, TopologyAttrs, ConnectivityAttrs
from MDAnalysis.core import selection as sel
from unittest.mock import patch

def test_TopologyParsers(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    assert top.lines == snapshot

def test_TopologyAttrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
        ta = TopologyAttrs(top.attrs)
    assert ta.lines == snapshot

def test_ConnectivityAttrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
        ca = ConnectivityAttrs(top.attrs)
    assert ca.lines == snapshot
