#!/usr/bin/env python
"""
Generate topology_defaults.txt:

A table of whether TopologyAttrs are atomwise, residuewise, or segmentwise, and their defaults
"""

from MDAnalysis.core.topologyattrs import AtomAttr, ResidueAttr, SegmentAttr 
from base import TOPOLOGY_CLS, write_rst_table

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

if __name__ == '__main__':
    lines = get_lines()
    write_rst_table(lines, HEADINGS, 'topology_defaults.txt')