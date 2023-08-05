#!/usr/bin/env python
"""
Generates:
    - ../formats/selection_exporters.txt
"""

from MDAnalysis import _SELECTION_WRITERS
from base import TableWriter
import helpers
from typing import Final

SELECTION_DESCRIPTIONS: Final = {
    'vmd': 'VMD macros, available in Representations',
    'pml': 'PyMOL selection string',
    'ndx': 'GROMACS index file',
    'str': 'CHARMM selection of individual atoms',
    'spt': 'Jmol selection commands',
}

def _program(klass) -> str:
    # classes have multiple formats.
    # First tends to be the program name, second is extension
    p = klass.format
    if isinstance(p, (list, tuple)):
        p = p[0]
    return helpers.sphinx_link(txt=p)

def _extension(klass) -> str:
    return klass.ext

def _description(klass) -> str:
    return SELECTION_DESCRIPTIONS[klass.ext]

def _class(klass) -> str:
    return helpers.sphinx_class(klass=klass, tilde=False)

class SelectionExporterWriter():
    def __init__(self):
        table_writer = TableWriter(
            headings = ['Program', 'Extension', 'Description', 'Class'],
            filename = 'formats/selection_exporter_formats.txt',
            include_table = 'Supported selection exporters',
            sort = True,
            inputs = set(_SELECTION_WRITERS.values()),
            methods = {
                'Program': _program,
                'Extension': _extension,
                'Description': _description,
                'Class': _class,
            }
        )
        table_writer.do_stuff()
        self.lines = table_writer.lines


if __name__ == '__main__':
    SelectionExporterWriter()
