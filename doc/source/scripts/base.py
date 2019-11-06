from __future__ import print_function
import os
import tabulate
from MDAnalysis import _TOPOLOGY_ATTRS

ignore = ('topologyattrs', 'atomattrs', 'residueattrs',
          'segmentattrs', 'indices', 'resindices', 'segindices')

TOPOLOGY_CLS = sorted(set([x for x in _TOPOLOGY_ATTRS.values()
                           if not x.attrname in ignore]), 
                      key=lambda x: x.attrname)

DESCRIPTIONS = {
    'CRD': 'CHARMM CARD file',
    'CONFIG': 'DL_Poly CONFIG file',
    'DCD': 'CHARMM, NAMD, or LAMMPS binary trajectory',
    'HISTORY': 'DL_Poly HISTORY file',
    'DMS': 'DESRES Molecular Structure file',
    'XPDB': 'Extended PDB file',
    'GMS': 'GAMESS file',
    'GRO': 'GROMACS structure file',
    'GSD': 'HOOMD GSD file',
    'XML': 'HOOMD XML file',
    'DATA': 'LAMMPS data file',
    'INPCRD': 'AMBER restart file',
    'LAMMPS': 'a LAMMPS DCD trajectory',
    'LAMMPSDUMP': 'LAMMPS ascii dump file',
    'MMTF' : 'MMTF file',
    'NCDF': 'AMBER NETCDF format',
    'MOL2': 'Tripos MOL2 file',
    'PDB': 'Standard PDB file',
    'PDBQT' : 'PDBQT file',
    'PQR' : 'PQR file',
    'PSF' : 'CHARMM, NAMD, or XPLOR PSF file',
    'TOP': 'AMBER topology file',
    'TPR': 'GROMACS run topology file',
    'TRJ': 'AMBER ASCII trajectories',
    'TRR': 'GROMACS TRR trajectory',
    'TRZ': 'IBIsCO or YASP binary trajectory',
    'TXYZ': 'Tinker file',
    'XTC': 'GROMACS compressed trajectory',
    'XYZ': 'XYZ file',
}

ATTR_DESCRIPTIONS = {
    'altLocs': 'Alternate location',
    'atomiccharges': 'Atomic number',
    'atomnums': '?',
    'bfactors': 'alias of tempfactor',
    'chainIDs': 'chain ID',
    'charges': 'partial atomic charge',
    'elements': 'atom element',
    'icodes': 'atom insertion code',
    'models': 'model number (from 0)',
    'molnums': '[molecules] number (from 0)',
    'moltypes': '[moleculetype] name',
    'names': 'atom names',
    'occupancies': 'atom occupancy',
    'radii': 'atomic radius',
    'record_types': 'ATOM / HETATM',
    'resnames': 'residue name (except GSD has ints)',
    'tempfactors': 'B-factor',
    'type_indices': 'amber atom type number',
}

ATTRS = {c.attrname:(c.singular, ATTR_DESCRIPTIONS.get(c.attrname, '')) for c in _TOPOLOGY_ATTRS.values()}

def sphinx_class(klass, tilde=True):
    prefix = '~' if tilde else ''
    return ':class:`{}{}.{}`'.format(prefix, klass.__module__, klass.__name__)

def sphinx_meth(meth, tilde=True):
    prefix = '~' if tilde else ''
    return ':meth:`{}{}.{}'.format(prefix, meth.__module__, meth.__qualname__)

def sphinx_ref(txt, label=None, add_label=True):
    if label is None:
        label = txt
    suffix = '-label' if add_label else ''
    return ':ref:`{} <{}{}>`'.format(txt, label, suffix)

def write_rst_table(lines=[], headings=[], filename='filename.txt'):
    path = os.getcwd()
    if not 'scripts' in path:
        filename = os.path.join('scripts', filename)
    with open(filename, 'w') as f:
        print(tabulate.tabulate(lines, headers=headings, tablefmt='rst'), file=f)
    print('Wrote ', filename)