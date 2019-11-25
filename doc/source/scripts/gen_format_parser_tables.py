#!/usr/bin/env python
"""
Generates:
    - ../formats/format_overview.txt : table of format overview
    - ../formats/coordinate_readers.txt: Coordinate reader information
    - ../formats/classes/*.txt : small tables of parser/reader/writer classes for each format
"""

from collections import defaultdict

from MDAnalysis import _READERS, _SINGLEFRAME_WRITERS, _PARSERS
from base import TableWriter
from core import DESCRIPTIONS

FILE_TYPES = defaultdict(dict)

for clstype, dct in (('Coordinate reader', _READERS),
                        ('Coordinate writer', _SINGLEFRAME_WRITERS),
                        ('Topology parser', _PARSERS)):
    for fmt, klass in dct.items():
        if fmt in ('CHAIN', 'MEMORY', 'MINIMAL', 'NULL'):
            continue  # get their own pages
        FILE_TYPES[fmt][clstype] = klass

sorted_types = sorted(FILE_TYPES.items())

class FormatOverview(TableWriter):
    filename = 'formats/format_overview.txt'
    preprocess = ['keys']
    headings = ['File type', 'Description', 'Topology', 'Coordinates', 'Read/Write']

    def _set_up_input(self):
        return sorted_types
    
    def _file_type(self, fmt, handlers):
        return self.sphinx_ref(fmt, self.keys[-1])
    
    def _keys(self, fmt, handlers):
        if fmt in DESCRIPTIONS:
            key = fmt
        else:
            key = list(handlers.values())[0].format[0]

            # raise an informative error
            if key not in DESCRIPTIONS:
                key = fmt
        return key
    
    def _description(self, fmt, handlers):
        return DESCRIPTIONS[self.keys[-1]]
    
    def _topology(self, fmt, handlers):
        if 'Topology parser' in handlers:
            return 'Yes'
        return ''
    
    def _coordinates(self, fmt, handlers):
        if 'Coordinate reader' in handlers:
            return 'Yes'
        return ''
    
    def _read_write(self, fmt, handlers):
        if 'Coordinate writer' in handlers:
            return 'Read and write'
        return 'Read only'

class CoordinateReaders(FormatOverview):
    filename = 'formats/coordinate_readers.txt'
    headings = ['File type', 'Description', 'Velocities', 'Forces']

    def _set_up_input(self):
        return [(x, y) for x, y in sorted_types if 'Coordinate reader' in y]

    def _velocities(self, fmt, handlers):
        if handlers['Coordinate reader'].units.get('velocity', None):
            return 'Yes'
        return ''
    
    def _forces(self, fmt, handlers):
        if handlers['Coordinate reader'].units.get('force', None):
            return 'Yes'
        return ''

class SphinxClasses(TableWriter):
    
    filename = 'formats/reference/classes/{}.txt'

    def __init__(self, fmt):
        self.filename = self.filename.format(fmt)
        self.fmt = fmt
        super(SphinxClasses, self).__init__()
    
    def get_lines(self):
        lines = []
        for label, klass in sorted(FILE_TYPES[self.fmt].items()):
            lines.append(['**{}**'.format(label), 
                          self.sphinx_class(klass, tilde=False)])
        self.lines = lines

if __name__ == '__main__':
    ov = FormatOverview()
    CoordinateReaders()
    for key in set(ov.fields['keys']):
        SphinxClasses(key)

        