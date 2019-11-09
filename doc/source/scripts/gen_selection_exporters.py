#!/usr/bin/env python
"""
Generates:
    - ../formats/selection_exporters.txt
"""

from MDAnalysis import _SELECTION_WRITERS
from base import TableWriter

SELECTION_DESCRIPTIONS = {
    'vmd': 'VMD macros, available in Representations',
    'pml': 'PyMOL selection string',
    'ndx': 'GROMACS index file',
    'str': 'CHARMM selection of individual atoms',
    'spt': 'Jmol selection commands',
}

class SelectionExporterWriter(TableWriter):
    headings = ['Program', 'Extension', 'Description', 'Class']
    filename = 'formats/selection_exporter_formats.txt'
    sort = True

    def _set_up_input(self):
        return set(_SELECTION_WRITERS.values())
    
    def _program(self, klass):
        # classes have multiple formats.
        # First tends to be the program name, second is extension
        p = klass.format
        if isinstance(p, (list, tuple)):
            p = p[0]
        return self.sphinx_link(p)
    
    def _extension(self, klass):
        return klass.ext
    
    def _description(self, klass):
        return SELECTION_DESCRIPTIONS[klass.ext]
    
    def _class(self, klass):
        return self.sphinx_class(klass, tilde=False)
    

if __name__ == '__main__':
    SelectionExporterWriter()
