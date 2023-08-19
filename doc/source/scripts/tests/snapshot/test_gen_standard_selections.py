from unittest.mock import patch

from gen_standard_selections import StandardSelectionTable
from MDAnalysis.core import selection as sel


def test_StandardSelectionTable_protein(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable(
            "protein",
            klass=sel.ProteinSelection,
            attribute_name="prot_res",
            sort=True,
        )
    assert ss.table_writer.lines == snapshot


def test_StandardSelectionTable_protein_backbone(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable(
            "protein_backbone",
            klass=sel.BackboneSelection,
            attribute_name="bb_atoms",
            sort=True,
        )
    assert ss.table_writer.lines == snapshot


def test_StandardSelectionTable_nucleic(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable(
            "nucleic",
            klass=sel.NucleicSelection,
            attribute_name="nucl_res",
            sort=True,
        )
    assert ss.table_writer.lines == snapshot


def test_StandardSelectionTable_nucleic_backbone(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable(
            "nucleic_backbone",
            klass=sel.NucleicBackboneSelection,
            attribute_name="bb_atoms",
            sort=True,
        )
    assert ss.table_writer.lines == snapshot


def test_StandardSelectionTable_base(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable(
            "base",
            klass=sel.BaseSelection,
            attribute_name="base_atoms",
            sort=True,
        )
    assert ss.table_writer.lines == snapshot


def test_StandardSelectionTable_nucleic_sugar(snapshot):
    with patch("builtins.open"):
        ss = StandardSelectionTable(
            "nucleic_sugar",
            klass=sel.NucleicSugarSelection,
            attribute_name="sug_atoms",
            sort=True,
        )
    assert ss.table_writer.lines == snapshot
