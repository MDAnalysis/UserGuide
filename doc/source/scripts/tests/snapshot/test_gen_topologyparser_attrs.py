from unittest.mock import patch

from gen_topologyparser_attrs import ConnectivityAttrs, TopologyAttrs, TopologyParsers


def test_TopologyParsers_lines(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    assert top.table_writer.lines == snapshot


def test_TopologyParsers_attrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    assert top.attrs == snapshot


def test_TopologyAttrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
        ta = TopologyAttrs(top.attrs)
    assert ta.table_writer.lines == snapshot


def test_ConnectivityAttrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
        ca = ConnectivityAttrs(top.attrs)
    assert ca.table_writer.lines == snapshot
