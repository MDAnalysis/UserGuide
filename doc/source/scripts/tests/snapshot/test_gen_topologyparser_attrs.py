from unittest.mock import patch

from gen_topologyparser_attrs import (
    ConnectivityAttrs,
    TopologyAttrs,
    TopologyParsers,
    get_format_attrs,
)


def test_TopologyParsers_lines(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    assert top.table_writer.lines == snapshot


def test_TopologyParsers_attrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
    attrs = get_format_attrs(top)
    assert attrs == snapshot


def test_TopologyAttrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
        attrs = get_format_attrs(top)
        ta = TopologyAttrs(attrs)
    assert ta.table_writer.lines == snapshot


def test_ConnectivityAttrs(snapshot):
    with patch("builtins.open"):
        top = TopologyParsers()
        attrs = get_format_attrs(top)
        ca = ConnectivityAttrs(attrs)
    assert ca.table_writer.lines == snapshot
