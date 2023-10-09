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
from typing import Any

from base import TableWriter
from MDAnalysis.core import selection as sel


def _chunk_list(lst: list[str], chunk_size: int = 8) -> list[list[str]]:
    return [lst[i : i + chunk_size] for i in range(0, len(lst), chunk_size)]


# override get_lines as there are no headings
def _custom_get_lines(
    *,
    klass: sel.Selection,
    attribute_name: str,
    sort: bool = False,
    chunk_size: int = 8
) -> list[list[str]]:
    selected = getattr(klass, attribute_name)
    if sort:
        selected = sorted(selected)

    table = _chunk_list(list(selected), chunk_size=chunk_size)
    return table


class StandardSelectionTable:
    def __init__(self, filename: str, *args: Any, **kwargs: Any) -> None:

        lines = _custom_get_lines(*args, **kwargs)
        self.table_writer = TableWriter(
            sort=False,
            filename="generated/selections/{}.txt".format(filename),
            lines=lines,
            column_spec=[],
        )
        self.table_writer.write_table()


if __name__ == "__main__":
    StandardSelectionTable(
        "protein",
        klass=sel.ProteinSelection,
        attribute_name="prot_res",
        sort=True,
    )
    StandardSelectionTable(
        "protein_backbone",
        klass=sel.BackboneSelection,
        attribute_name="bb_atoms",
    )
    StandardSelectionTable(
        "nucleic", klass=sel.NucleicSelection, attribute_name="nucl_res"
    )
    StandardSelectionTable(
        "nucleic_backbone",
        klass=sel.NucleicBackboneSelection,
        attribute_name="bb_atoms",
    )
    StandardSelectionTable(
        "base", klass=sel.BaseSelection, attribute_name="base_atoms"
    )
    StandardSelectionTable(
        "nucleic_sugar",
        klass=sel.NucleicSugarSelection,
        attribute_name="sug_atoms",
    )
