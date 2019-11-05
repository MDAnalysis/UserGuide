#!/usr/bin/env python
"""
Generate topology_defaults.txt:

A table of whether TopologyAttrs are atomwise, residuewise, or segmentwise, and their defaults
"""

from __future__ import print_function
import os
import sys
import tabulate
from collections import defaultdict
from MDAnalysis import _TOPOLOGY_ATTRS
from MDAnalysis.core.topologyattrs import AtomAttr, ResidueAttr, SegmentAttr 

ignore = ('topologyattrs', 'atomattrs', 'residueattrs', 'segmentattrs', 'indices', 'resindices', 'segindices')

TOPOLOGY_CLS = sorted(set([x for x in _TOPOLOGY_ATTRS.values()
                           if not x.attrname in ignore]), 
                      key=lambda x: x.attrname)

HEADINGS = ('Atom', 'AtomGroup', 'default', 'level', 'type')

DEFAULTS = {
    'resids': 'continuous sequence from 1 to n_residues',
    'resnums': 'continuous sequence from 1 to n_residues',
    'ids': 'continuous sequence from 1 to n_atoms',
}

for klass in TOPOLOGY_CLS:
    if klass.attrname not in DEFAULTS:
        try:
            DEFAULTS[klass.attrname] = repr(klass._gen_initial_values(1, 1, 1)[0])
        except NotImplementedError:
            DEFAULTS[klass.attrname] = 'No default values'

def get_lines():
    lines = []
    for klass in TOPOLOGY_CLS:
        if issubclass(klass, AtomAttr):
            level = 'atom'
        elif issubclass(klass, ResidueAttr):
            level = 'residue'
        elif issubclass(klass, SegmentAttr):
            level = 'segment'
        lines.append((klass.attrname, klass.singular, DEFAULTS[klass.attrname], level, klass.dtype))
    return lines

def write_table(lines):
    filename = 'topology_defaults.txt'
    path = os.getcwd()
    if not 'scripts' in path:
        filename = os.path.join('scripts', filename)
    with open(filename, 'w') as f:
        print(tabulate.tabulate(lines, headers=HEADINGS, tablefmt='rst'), file=f)

if __name__ == '__main__':
    lines = get_lines()
    write_table(lines)