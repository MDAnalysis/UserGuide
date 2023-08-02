
from doc.source.scripts.gen_standard_selections import StandardSelectionTable
from MDAnalysis.core import selection as sel


def test_StandardSelectionTable_protein(snapshot):
    ss = StandardSelectionTable('protein', sel.ProteinSelection, 'prot_res', True)
    assert ss.lines == snapshot


def test_StandardSelectionTable_protein_backbone(snapshot):
    ss = StandardSelectionTable('protein_backbone', sel.BackboneSelection, 'bb_atoms', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_nucleic(snapshot):
    ss = StandardSelectionTable('nucleic', sel.NucleicSelection, 'nucl_res', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_nucleic_backbone(snapshot):
    ss = StandardSelectionTable('nucleic_backbone', sel.NucleicBackboneSelection, 'bb_atoms', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_base(snapshot):
    ss = StandardSelectionTable('base', sel.BaseSelection, 'base_atoms', True)
    assert ss.lines == snapshot

def test_StandardSelectionTable_nucleic_sugar(snapshot):
    ss = StandardSelectionTable('nucleic_sugar', sel.NucleicSugarSelection, 'sug_atoms', True)
    assert ss.lines == snapshot
