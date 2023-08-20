from unittest.mock import patch

from gen_topologyparser_attrs import ConnectivityAttrs, TopologyAttrs, TopologyParsers
from MDAnalysis.core import selection as sel


def test_TopologyParsers_lines(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    assert top.lines == snapshot


def test_TopologyParsers_attrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    assert top.attrs == snapshot


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
