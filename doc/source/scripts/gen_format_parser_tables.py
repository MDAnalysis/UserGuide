#!/usr/bin/env python
"""
Generates:
    - ../formats/format_overview.txt : table of format overview
    - ../formats/coordinate_readers.txt: Coordinate reader information
    - ../formats/classes/*.txt : small tables of parser/reader/writer classes for each format
"""

from __future__ import print_function
import os
import tabulate
from collections import defaultdict

from MDAnalysis import _READERS, _SINGLEFRAME_WRITERS, _PARSERS
from base import write_rst_table, sphinx_class, sphinx_ref, DESCRIPTIONS

overview_headings = ['File type', 'Description', 'Topology', 'Coordinates', 'Read/Write']
coordinate_headings = ['File type', 'Velocities', 'Forces']
row_labels = ['Topology parser', 'Coordinate reader', 'Coordinate writer']


FILE_TYPES = defaultdict(dict)

for clstype, dct in (('Coordinate reader', _READERS),
                        ('Coordinate writer', _SINGLEFRAME_WRITERS),
                        ('Topology parser', _PARSERS)):
    for fmt, klass in dct.items():
        if fmt in ('CHAIN', 'MEMORY', 'MINIMAL', 'NULL'):
            continue  # get their own pages
        FILE_TYPES[fmt][clstype] = klass


def get_overview():
    overview = []
    coordinates = []
    non_repeats = set()
    for fmt, handlers in sorted(FILE_TYPES.items()):
        try:
            desc = DESCRIPTIONS[fmt]
            key = fmt
        except KeyError:
            key = list(handlers.values())[0].format[0]
            desc = DESCRIPTIONS[key]

        non_repeats.add(key)
        top = ('', 'Yes')[int('Topology parser' in handlers)]
        coords = ('', 'Yes')[int('Coordinate reader' in handlers)]
        rw = ('Read only', 'Read and write')[int('Coordinate writer' in handlers)]
        ref = sphinx_ref(fmt, key)
        overview.append((ref, desc, top, coords, rw))

        try:
            reader = handlers['Coordinate reader']
        except KeyError:
            continue
        v = 'Yes' if reader.units.get('velocity', None) else ''
        f = 'Yes' if reader.units.get('force', None) else ''
        coordinates.append((ref, v, f))
    return overview, coordinates, non_repeats


def get_classes(fmts):
    classes = defaultdict(list)
    for fmt in fmts:
        handlers = FILE_TYPES[fmt]
        for label, klass in sorted(handlers.items()):
            classes[fmt].append(['**{}**'.format(label), sphinx_class(klass, tilde=False)])
    return classes


if __name__ == '__main__':
    overview, coordinates, non_repeats = get_overview()
    write_rst_table(lines=overview, headings=overview_headings,
                    filename='../formats/format_overview.txt')

    write_rst_table(lines=coordinates, headings=coordinate_headings,
                    filename='../formats/coordinate_readers.txt')
    classes = get_classes(non_repeats)
    for fmt, lines in classes.items():
        write_rst_table(lines, headings=[],
                        filename='../formats/classes/{}.txt'.format(fmt))


        