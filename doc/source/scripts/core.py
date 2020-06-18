from __future__ import print_function
import os
import tabulate
from MDAnalysis import _TOPOLOGY_ATTRS

# ====== TOPOLOGY ====== #

DESCRIPTIONS = {
    'CHEMFILES': 'chemfiles object',
    'PARMED': 'parmed Structure',
    'CRD': 'CHARMM CARD file',
    'CONFIG': 'DL_Poly CONFIG file',
    'COOR': 'NAMDBIN file',
    'DCD': 'CHARMM, NAMD, or LAMMPS binary trajectory',
    'HISTORY': 'DL_Poly HISTORY file',
    'DMS': 'DESRES Molecular Structure file',
    'XPDB': 'Extended PDB file',
    'FHIAIMS': 'FHIAIMS file',
    'GMS': 'GAMESS file',
    'GRO': 'GROMACS structure file',
    'GSD': 'HOOMD GSD file',
    'XML': 'HOOMD XML file',
    'DATA': 'LAMMPS data file',
    'INPCRD': 'AMBER restart file',
    'ITP': 'GROMACS portable topology file',
    'IN': 'FHIAIMS file',
    'LAMMPS': 'a LAMMPS DCD trajectory',
    'LAMMPSDUMP': 'LAMMPS ascii dump file',
    'MMTF' : 'MMTF file',
    'NAMDBIN': 'NAMDBIN file',
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

base_attrnames = set(['atomattrs', 'residueattrs', 'segmentattrs', 'topologyattrs'])

core_attrnames = set(['indices', 'resindices', 'segindices'])

BASE_ATTRS = {k:v for k, v in ATTRS.items() if k in base_attrnames}

NON_BASE_ATTRS = {k:v for k, v in ATTRS.items() if k not in base_attrnames}

NON_CORE_ATTRS = {k:v for k, v in NON_BASE_ATTRS.items() if k not in core_attrnames}

TOPOLOGY_CLS = sorted(set([x for x in _TOPOLOGY_ATTRS.values()
                           if x.attrname in NON_CORE_ATTRS.keys()]), 
                      key=lambda x: x.attrname)