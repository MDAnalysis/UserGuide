#!/usr/bin/env python
"""
Generate topology_defaults.txt:

A table of whether TopologyAttrs are atomwise, residuewise, or segmentwise, and their defaults
"""

from MDAnalysis.core.topologyattrs import AtomAttr, ResidueAttr, SegmentAttr
from core import TOPOLOGY_CLS
from base import TableWriter


DEFAULTS = {
    'resids': 'continuous sequence from 1 to n_residues',
    'resnums': 'continuous sequence from 1 to n_residues',
    'ids': 'continuous sequence from 1 to n_atoms',
}

class TopologyDefaults(TableWriter):
    filename = 'generated/topology/defaults.txt'
    headings = ('Atom', 'AtomGroup', 'default', 'level', 'type')
    sort = True

    def _set_up_input(self):
        return TOPOLOGY_CLS

    def _atom(self, klass):
        return klass.attrname

    def _atomgroup(self, klass):
        return klass.singular

    def _default(self, klass):
        try:
            return DEFAULTS[klass.attrname]
        except KeyError:
            try:
                return repr(klass._gen_initial_values(1, 1, 1)[0])
            except NotImplementedError:
                return 'No default values'

    def _level(self, klass):
        if issubclass(klass, AtomAttr):
            level = 'atom'
        elif issubclass(klass, ResidueAttr):
            level = 'residue'
        elif issubclass(klass, SegmentAttr):
            level = 'segment'
        else:
            raise ValueError
        return level

    def _type(self, klass):
        return klass.dtype

if __name__ == '__main__':
    TopologyDefaults()
