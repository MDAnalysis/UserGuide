#!/usr/bin/env python
"""
Writes standard selection names for:
    - protein residues
    - protein backbone atoms
    - nucleic residues
    - nucleic backbone atoms
    - nucleobase atoms
    - nucleic sugar atoms
"""
from base import TableWriter
from MDAnalysis.core import selection as sel


def chunk_list(lst, n=8):
    return [lst[i : i + n] for i in range(0, len(lst), n)]


class StandardSelectionTable(TableWriter):
    sort = False
    filename = "generated/selections/{}.txt"

    def __init__(self, filename, *args, **kwargs):
        self.filename = self.filename.format(filename)
        super(StandardSelectionTable, self).__init__(*args, **kwargs)

    # override get_lines as there are no headings
    def get_lines(self, klass, attr, sort=False, n=8):
        selected = getattr(klass, attr)
        if sort:
            selected = sorted(selected)

        table = chunk_list(list(selected), n=n)
        self.lines = table


if __name__ == "__main__":
    StandardSelectionTable("protein", sel.ProteinSelection, "prot_res", True)
    StandardSelectionTable(
        "protein_backbone", sel.BackboneSelection, "bb_atoms"
    )
    StandardSelectionTable("nucleic", sel.NucleicSelection, "nucl_res")
    StandardSelectionTable(
        "nucleic_backbone", sel.NucleicBackboneSelection, "bb_atoms"
    )
    StandardSelectionTable("base", sel.BaseSelection, "base_atoms")
    StandardSelectionTable(
        "nucleic_sugar", sel.NucleicSugarSelection, "sug_atoms"
    )
