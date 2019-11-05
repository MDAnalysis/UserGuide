#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import tabulate
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

DESCRIPTIONS = {
    'CRD': 'CHARMM CARD file',
    'CONFIG': 'DL_Poly CONFIG file',
    'HISTORY': 'DL_Poly HISTORY file',
    'DMS': 'DESRES Molecular Structure file',
    'XPDB': 'Extended PDB file',
    'GMS': 'GAMESS file',
    'GRO': 'GROMACS structure file',
    'GSD': 'HOOMD GSD file',
    'XML': 'HOOMD XML file',
    'DATA': 'LAMMPS data file',
    'LAMMPSDUMP': 'LAMMPS ascii dump file',
    'MMTF' : 'MMTF file',
    'MOL2': 'Tripos MOL2 file',
    'PDB': 'Standard PDB file',
    'PDBQT' : 'PDBQT file',
    'PQR' : 'PQR file',
    'PSF' : 'CHARMM, NAMD, or XPLOR PSF file',
    'TOP': 'AMBER topology file',
    'TPR': 'GROMACS run topology file',
    'TXYZ': 'Tinker file',
    'XYZ': 'XYZ file',
}

COORDINATES = ('CRD', 'CONFIG', 'HISTORY', 'DMS', 
               'GMS', 'GRO', 'GSD', 'DATA', 
               'LAMMPSDUMP', 'MMTF', 'MOL2', 
               'PDB', 'PDBQT', 'PQR', 'TXYZ', 'XYZ')

NOTES =  {
    'TPR': 'indexes atom ids and resids from 0',
    'DMS': 'indexes atom ids from 0',
    'PSF': 'atom ids are off by 1 from the file; MDAnalysis subtracts 1 from the file atom ids',

}

MANDATORY_ATTRS = set(mandatory_attrs)

HEADINGS = ('Format', 'Type',
            'Attributes read', 'Attributes guessed',
            # 'Contains coordinates',  # too long
            # 'Notes',  # too long
            'Class')

def get_lines():
    lines = []
    for p in PARSER_TESTS:
        f = p.parser.format
        if isinstance(f, (list, tuple)):
            key = f[0]
            fstr = ':ref:`{} <{}>`'.format(', '.join(f), '{}-label'.format(key))
        else:
            key = f
            fstr = ':ref:`{} <{}>`'.format(f, '{}-label'.format(f))
        
        
        desc = DESCRIPTIONS[key]  # description
        read = sorted(set(p.expected_attrs) - MANDATORY_ATTRS)  # attributes read
        guessed = sorted(p.guessed_attrs)  # attributes guessed
        # coordinates = key in COORDINATES  # contains coordinates
        # notes = NOTES.get(key, '')  # other notes

        if key in COORDINATES:
            fstr += '\ [#coordinates]_'
        klass = ':class:`~'+'.'.join([p.parser.__module__, p.parser.__name__])+'`'  # class

        lines.append((fstr, desc, ', '.join(read), ', '.join(guessed), klass))
    return lines

def write_topology_table():
    lines = get_lines()
    filename = 'topology_parsers.txt'
    path = os.getcwd()
    if not 'scripts' in path:
        filename = os.path.join('scripts', filename)
    with open(filename, 'w') as f:
        print(tabulate.tabulate(lines, headers=HEADINGS, tablefmt='rst'), file=f)

if __name__ == '__main__':
    write_topology_table()