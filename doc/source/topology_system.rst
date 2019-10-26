.. -*- coding: utf-8 -*-
.. _topology-system:

=====================
The topology system
=====================

MDAnalysis groups static data about a :class:`~MDAnalysis.core.universe.Universe` into its topology. This is typically loaded from a topology file. Topology information falls into 3 categories:

    * :ref:`Atom containers (Residues and Segments) <residues-and-segments>`
    * :ref:`Atom attributes (e.g. name, mass, bfactor) <topology-attributes>`
    * :ref:`Topology objects: bonds, angles, dihedrals, impropers <topology-objects>`

Users will almost never interact directly with a :class:`~MDAnalysis.core.topology.Topology`. Modifying atom containers or topology attributes is typically done through :class:`~MDAnalysis.core.universe.Universe`. Methods for viewing containers or topology attributes, or for calculating topology object values, are accessed through :class:`~MDAnalysis.core.groups.AtomGroup`.


.. _topology-attributes:

Topology attributes
===================

MDAnalysis supports a range of topology attributes for each :class:`~MDAnalysis.core.groups.Atom` and :class:`~MDAnalysis.core.groups.AtomGroup`. If an attribute is defined for an Atom, it will be for an AtomGroup, and vice versa -- however, they are accessed with singular and plural versions of the attribute specifically.

---------------------------------
Canonical attributes
---------------------------------

These attributes are derived for every :class:`~MDAnalysis.core.universe.Universe`, including Universes created with :meth:`~MDAnalysis.core.universe.Universe.empty`. They encode the MDAnalysis order of each object.

+----------+---------------+---------------------------------------------+
| **Atom** | **AtomGroup** | **Description**                             |
+----------+---------------+---------------------------------------------+
| index    | indices       | MDAnalysis canonical atom index (from 0)    |
+----------+---------------+---------------------------------------------+
| resindex | resindices    | MDAnalysis canonical residue index (from 0) |
+----------+---------------+---------------------------------------------+
| segindex | segindices    | MDAnalysis segment index (from 0)           |
+----------+---------------+---------------------------------------------+

The following attributes are read or guessed from every format supported by MDAnalysis.

+----------+---------------+---------------------------------------------------+
| **Atom** | **AtomGroup** | **Description**                                   |
+----------+---------------+---------------------------------------------------+
| id       | ids           | atom serial (from 1, except PSF/DMS/TPR formats)  |
+----------+---------------+---------------------------------------------------+
| mass     | masses        | atom mass (guessed, default: 0.0)                 |
+----------+---------------+---------------------------------------------------+
| resid    | resids        | residue number (from 1, except for TPR)           |
+----------+---------------+---------------------------------------------------+
| resnum   | resnums       | alias of resid                                    |
+----------+---------------+---------------------------------------------------+
| segid    | segids        | names of segments (default: 'SYSTEM')             |
+----------+---------------+---------------------------------------------------+
| type     | types         | atom name, atom element, or force field atom type |
+----------+---------------+---------------------------------------------------+

---------------------------------
Format-specific attributes
---------------------------------

+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| **Atom**     | **AtomGroup** | **Description**                               | **Supported formats**                                                                                     |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| altLoc       | altLocs       | Alternate location                            | MMTF, PDB, PDBQT, PDBx                                                                                    |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| atomiccharge | atomiccharges | Atomic number                                 | GAMESS                                                                                                    |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| atomnum      | atomnums      | ?                                             | DMS                                                                                                       |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| bfactor      | bfactors      | alias of tempfactor                           | MMTF                                                                                                      |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| chainID      | chainIDs      | chain ID                                      | DMS                                                                                                       |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| charge       | charges       | partial atomic charge                         | DMS, HOOMD GSD, HOOMD XML, LAMMPS, MMTF, MOL2, PDBQT, PSF, PRMTOP, TPR                                    |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| element      | elements      | atom element                                  | PRMTOP                                                                                                    |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| icode        | icodes        | atom insertion code                           | MMTF, PDB, PQR                                                                                            |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| model        | models        | model number (from 0)                         | MMTF                                                                                                      |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| molnum       | molnums       | [ molecules ] number (from 0)                 | TPR                                                                                                       |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| moltype      | moltypes      | [ moleculetype ] name                         | TPR                                                                                                       |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| name         | names         | atom names                                    | CRD, DL Poly, DMS, GAMESS, GRO, HOOMD GSD, MMTF, MOL2, PDB, PDBQT, PDBx, PQR, PSF, PRMTOP, TPR, TXYZ, XYZ |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| occupancy    | occupancies   | atom occupancy                                | MMTF, PDBx, PDB, PDBQT                                                                                    |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| radius       | radii         | atomic radius                                 | GSD, HOOMD XML, PQR                                                                                       |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| record_type  | record_types  | ATOM / HETATM                                 | PDB, PDBQT, PQR, PDBx                                                                                     |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| resname      | resnames      | residue name (except GSD, which has integers) | CRD, DMS, GRO, HOOMD GSD, MMTF, MOL2, PDB, PDBQT, PDBx, PQR, PSF, PRMTOP, TPR                             |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| tempfactor   | tempfactors   | B-factor                                      | CRD, PDB, PDBx, PDBQT                                                                                     |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+
| type_index   | type_indices  | amber atom type number                        | PRMTOP                                                                                                    |
+--------------+---------------+-----------------------------------------------+-----------------------------------------------------------------------------------------------------------+

---------------------------------
Adding or modifying TopologyAttrs
---------------------------------


.. _topology-objects:

Topology objects
================

MDAnalysis defines four types of :class:`~MDAnalysis.core.topologyobjects.TopologyObject` by connectivity:

    * bonds
    * angles
    * dihedrals
    * improper angles

----------------------------
Generating
----------------------------

TopologyObjects can only be read from these file formats:

    * PSF
    * GSD
    * TPR
    * PDB (with CONECT records)
    * Prmtop
    * DMS
    * mol2
    * LAMMPS data
    * Tinker TXYZ
    * Hoomd XML
    * MMTF

Bonds can be guessed based on distance and Van der Waals' radii with :meth:`AtomGroup.guess_bonds <MDAnalysis.core.groups.AtomGroup.guess_bonds>`.


Users can also define new TopologyObjects through an :class:`~MDAnalysis.core.groups.AtomGroup`.



----------------------------
Calculating values
----------------------------

