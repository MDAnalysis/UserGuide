#!/usr/bin/env python
"""
Generate three tables:
    - connectivityattrs.txt: A table of supported formats for bonds, angles, dihedrals, impropers
    - topologyattrs.txt: A table of supported formats for non-connectivity attributes.
    - topology_parsers.txt: A table of all formats and the attributes they read and guess

This script imports the testsuite, which tests these.
"""
import os
import sys
from collections import defaultdict
from core import DESCRIPTIONS, NON_CORE_ATTRS
from base import TableWriter

from MDAnalysisTests.topology.base import mandatory_attrs
from MDAnalysisTests.topology.test_crd import TestCRDParser
from MDAnalysisTests.topology.test_dlpoly import TestDLPHistoryParser, TestDLPConfigParser
from MDAnalysisTests.topology.test_dms import TestDMSParser
from MDAnalysisTests.topology.test_fhiaims import TestFHIAIMS
from MDAnalysisTests.topology.test_gms import GMSBase
from MDAnalysisTests.topology.test_gro import TestGROParser
from MDAnalysisTests.topology.test_gsd import TestGSDParser
from MDAnalysisTests.topology.test_hoomdxml import TestHoomdXMLParser
from MDAnalysisTests.topology.test_lammpsdata import LammpsBase, TestDumpParser
from MDAnalysisTests.topology.test_mmtf import TestMMTFParser
from MDAnalysisTests.topology.test_mol2 import TestMOL2Base
from MDAnalysisTests.topology.test_pdb import TestPDBParser
from MDAnalysisTests.topology.test_pdbqt import TestPDBQT
from MDAnalysisTests.topology.test_pqr import TestPQRParser
from MDAnalysisTests.topology.test_psf import PSFBase
from MDAnalysisTests.topology.test_top import TestPRMParser
from MDAnalysisTests.topology.test_tprparser import TPRAttrs
from MDAnalysisTests.topology.test_txyz import TestTXYZParser
from MDAnalysisTests.topology.test_xpdb import TestXPDBParser
from MDAnalysisTests.topology.test_xyz import XYZBase


PARSER_TESTS = (TestCRDParser, TestDLPHistoryParser, TestDLPConfigParser,
                TestDMSParser, TestFHIAIMS, GMSBase,
                TestGROParser, TestGSDParser, TestHoomdXMLParser,
                LammpsBase, TestMMTFParser, TestMOL2Base,
                TestPDBParser, TestPDBQT, TestPQRParser, PSFBase,
                TestPRMParser, TPRAttrs, TestTXYZParser,
                TestXPDBParser, XYZBase, TestDumpParser)


MANDATORY_ATTRS = set(mandatory_attrs)


parser_attrs = {}


for p in PARSER_TESTS:
    e, g = set(p.expected_attrs)-MANDATORY_ATTRS, set(p.guessed_attrs)
    # clunky hack for PDB
    if p is TestPDBParser:
        e.add('elements')
    parser_attrs[p.parser] = (e, g)


class TopologyParsers(TableWriter):
    headings = ['Format', 'Description', 'Attributes read', 'Attributes guessed']
    preprocess = ['keys']
    filename = 'formats/topology_parsers.txt'
    include_table = 'Table of supported topology parsers and the attributes read'
    sort = True

    def __init__(self):
        self.attrs = defaultdict(set)
        super(TopologyParsers, self).__init__()

    def _set_up_input(self):
        return [[x, *y] for x, y in parser_attrs.items()]

    def get_line(self, parser, expected, guessed):
        line = super(TopologyParsers, self).get_line(parser, expected, guessed)
        for a in expected|guessed:
            self.attrs[a].add(self.fields['Format'][-1])
        return line

    def _keys(self, parser, *args):
        f = parser.format
        if isinstance(f, (list, tuple)):
            key = f[0]
            label = ', '.join(f)
        else:
            key = label = f
        return (key, label)


    def _description(self, *args):
        key, label = self.keys[-1]
        return DESCRIPTIONS[key]


    def _format(self, *args):
        key, label = self.keys[-1]
        return self.sphinx_ref(label, key, suffix='-format')

    def _attributes_read(self, parser, expected, guessed):
        vals = sorted(expected - guessed)
        return ', '.join(vals)

    def _attributes_guessed(self, parser, expected, guessed):
        return ', '.join(sorted(guessed))


class TopologyAttrs(TableWriter):

    headings = ('Atom', 'AtomGroup', 'Description', 'Supported formats')
    filename = 'generated/topology/topologyattrs.txt'

    def __init__(self, attrs):
        self.attrs = attrs
        super(TopologyAttrs, self).__init__()

    def _set_up_input(self):
        return sorted([x, *y] for x, y in NON_CORE_ATTRS.items() if x not in MANDATORY_ATTRS)

    def _atom(self, name, singular, *args):
        return singular

    def _atomgroup(self, name, *args):
        return name

    def _description(self, name, singular, description):
        return description

    def _supported_formats(self, name, singular, description):
        return ', '.join(sorted(self.attrs[name]))


class ConnectivityAttrs(TopologyAttrs):
    headings = ('Atom', 'AtomGroup', 'Supported formats')
    filename = 'generated/topology/connectivityattrs.txt'

    def _set_up_input(self):
        inp = [[x]*3 for x in 'bonds angles dihedrals impropers'.split()]
        return inp


if __name__ == '__main__':
    top = TopologyParsers()
    TopologyAttrs(top.attrs)
    ConnectivityAttrs(top.attrs)
