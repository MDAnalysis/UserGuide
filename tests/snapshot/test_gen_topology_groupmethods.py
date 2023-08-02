
from doc.source.scripts.gen_topology_groupmethods import TransplantedMethods
from MDAnalysis.core import selection as sel


def test_TransplantedMethods(snapshot):
    tm = TransplantedMethods()
    assert tm.lines == snapshot
