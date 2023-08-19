#!/usr/bin/env python
"""
Generates:
    - ../formats/selection_exporters.txt
"""

import base
from base import TableWriter
from MDAnalysis import _SELECTION_WRITERS
from MDAnalysis.selections.base import SelectionWriterBase

SELECTION_DESCRIPTIONS = {
    "vmd": "VMD macros, available in Representations",
    "pml": "PyMOL selection string",
    "ndx": "GROMACS index file",
    "str": "CHARMM selection of individual atoms",
    "spt": "Jmol selection commands",
}


class SelectionExporterWriter:
    def __init__(self) -> None:
        def _program(klass: SelectionWriterBase) -> str:
            # classes have multiple formats.
            # First tends to be the program name, second is extension
            p = klass.format
            if isinstance(p, (list, tuple)):
                p = p[0]
            return base.sphinx_link(txt=p)

        def _extension(klass: SelectionWriterBase) -> str:
            return klass.ext  # type: ignore

        def _description(klass: SelectionWriterBase) -> str:
            return SELECTION_DESCRIPTIONS[klass.ext]

        def _class(klass: SelectionWriterBase) -> str:
            return base.sphinx_class(klass=klass, tilde=False)

        self.table_writer = TableWriter(
            filename="formats/selection_exporter_formats.txt",
            include_table="Supported selection exporters",
            sort=True,
            input_items=set(_SELECTION_WRITERS.values()),
            column_spec=[
                ("Program", _program),
                ("Extension", _extension),
                ("Description", _description),
                ("Class", _class),
            ],
            lines=[],
        )
        self.table_writer.generate_lines_and_write_table()


if __name__ == "__main__":
    SelectionExporterWriter()
