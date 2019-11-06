#!/usr/bin/env python
"""
Generate two tables:
    - connectivityattrs.txt: A table of supported formats for bonds, angles, dihedrals, impropers
    - topologyattrs.txt: A table of supported formats for non-connectivity attributes.
    - topology_parsers.txt: A table of all formats and the attributes they read and guess

This script imports the testsuite, which tests these.
"""
from __future__ import print_function
import os
import sys
import tabulate
from collections import defaultdict
from base import DESCRIPTIONS, write_rst_table, sphinx_class, sphinx_ref, ATTRS

import MDAnalysis as mda

tests = os.path.join(mda.__path__[0], '..', '..', 'testsuite', 'MDAnalysisTests')
sys.path.append(tests)

from topology.base import mandatory_attrs
from topology.test_crd import TestCRDParser
from topology.test_dlpoly import TestDLPHistoryParser, TestDLPConfigParser
from topology.test_dms import TestDMSParser
from topology.test_gms import GMSBase
from topology.test_gro import TestGROParser
from topology.test_gsd import TestGSDParser
from topology.test_hoomdxml import TestHoomdXMLParser
from topology.test_lammpsdata import LammpsBase, TestDumpParser
from topology.test_mmtf import TestMMTFParser
from topology.test_mol2 import TestMOL2Base
from topology.test_pdb import TestPDBParser
from topology.test_pdbqt import TestPDBQT
from topology.test_pqr import TestPQRParser
from topology.test_psf import PSFBase
from topology.test_top import TestPRMParser
from topology.test_tprparser import TPRAttrs
from topology.test_txyz import TestTXYZParser
from topology.test_xpdb import TestXPDBParser
from topology.test_xyz import XYZBase

PARSER_TESTS = (TestCRDParser, TestDLPHistoryParser, TestDLPConfigParser, 
                TestDMSParser, GMSBase,
                TestGROParser, TestGSDParser, TestHoomdXMLParser, 
                LammpsBase, TestMMTFParser, TestMOL2Base, 
                TestPDBParser, TestPDBQT, TestPQRParser, PSFBase, 
                TestPRMParser, TPRAttrs, TestTXYZParser, 
                TestXPDBParser, XYZBase, TestDumpParser)

MANDATORY_ATTRS = set(mandatory_attrs)


def get_lines():
    lines = []
    all_attrs = defaultdict(list)
    for p in PARSER_TESTS:
        f = p.parser.format
        if isinstance(f, (list, tuple)):
            key = f[0]
            fstr = sphinx_ref(', '.join(f), key)
        else:
            key = f
            fstr = sphinx_ref(f)        
        
        # desc = DESCRIPTIONS[key]  # description
        read = sorted((set(p.expected_attrs)-MANDATORY_ATTRS)-set(p.guessed_attrs))  # attributes read
        guessed = sorted(p.guessed_attrs)  # attributes guessed

        for a in set(read+guessed):
            all_attrs[a].append(fstr)

        lines.append((fstr, 
                    #   desc, 
                      ', '.join(read), 
                      ', '.join(guessed)))
    return lines, all_attrs

FORMAT_HEADINGS = ('Format', 
                #    'Type', 
                   'Attributes read', 
                   'Attributes guessed')
ATTR_HEADINGS = ('Atom', 'AtomGroup', 'Description', 'Supported formats')
CONNECTIVITY_HEADINGS = ('Atom', 'AtomGroup', 'Supported formats')

def get_attr_lines(attrs):
    lines = []
    for a, classes in sorted(attrs.items()):
        try:
            single, desc = ATTRS[a]
            lines.append((single, a, desc, ', '.join(sorted(classes))))
        except KeyError:
            pass

    bond_lines = []
    for a in ('bonds', 'angles', 'dihedrals', 'impropers'):
        classes = attrs.get(a, [])
        bond_lines.append((a, a, ', '.join(sorted(classes))))
    
    return lines, bond_lines


if __name__ == '__main__':
    format_lines, attrs = get_lines()
    write_rst_table(format_lines, headings=FORMAT_HEADINGS, filename='../formats/topology_parsers.txt')
    attr_lines, bond_lines = get_attr_lines(attrs)
    write_rst_table(attr_lines, headings=ATTR_HEADINGS, filename='topologyattrs.txt')
    write_rst_table(bond_lines, headings=CONNECTIVITY_HEADINGS, filename='connectivityattrs.txt')