from unittest.mock import patch

from gen_topology_groupmethods import TransplantedMethods
from MDAnalysis.core import selection as sel


def test_TransplantedMethods(snapshot):
    with patch("builtins.open"):
        tm = TransplantedMethods()
    assert tm.lines == snapshot
