
from doc.source.scripts.gen_standard_selections import StandardSelectionTable
from MDAnalysis.core import selection as sel
from unittest.mock import patch

def test_StandardSelectionTable_protein(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable('protein', sel.ProteinSelection, 'prot_res', True)
    assert ss.lines == snapshot


def test_StandardSelectionTable_protein_backbone(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable('protein_backbone', sel.BackboneSelection, 'bb_atoms', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_nucleic(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable('nucleic', sel.NucleicSelection, 'nucl_res', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_nucleic_backbone(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable('nucleic_backbone', sel.NucleicBackboneSelection, 'bb_atoms', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_base(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable('base', sel.BaseSelection, 'base_atoms', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_nucleic_sugar(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable('nucleic_sugar', sel.NucleicSugarSelection, 'sug_atoms', True)
    assert ss.lines == snapshot
