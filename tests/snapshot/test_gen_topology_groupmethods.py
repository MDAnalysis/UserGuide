
from doc.source.scripts.gen_topology_groupmethods import TransplantedMethods
from MDAnalysis.core import selection as sel
from unittest.mock import patch

def test_TransplantedMethods(snapshot):
    with patch("builtins.open"):
        tm = TransplantedMethods()
    assert tm.lines == snapshot
