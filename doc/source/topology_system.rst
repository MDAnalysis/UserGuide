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

--------------------
Canonical attributes
--------------------

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

.. _topologyobject-attr-label:

-----------------------------------
Bonds, angles, dihedrals, impropers
-----------------------------------

+-----------+---------------+-------------------------------------------------------------------------------------------------------+
| **Atom**  | **AtomGroup** | **Supported formats**                                                                                 |
+-----------+---------------+-------------------------------------------------------------------------------------------------------+
| bonds     | bonds         | DMS, GSD, HOOMD XML, LAMMPS data, MMTF, MOL2, PDB (with CONECT records), PRMTOP, PSF, Tinker XYZ, TPR |
+-----------+---------------+-------------------------------------------------------------------------------------------------------+
| angles    | angles        | GSD, HOOMD XML, LAMMPS data, PRMTOP, PSF, TPR                                                         |
+-----------+---------------+-------------------------------------------------------------------------------------------------------+
| dihedrals | dihedrals     | GSD, HOOMD XML, LAMMPS data, PRMTOP, PSF, TPR                                                         |
+-----------+---------------+-------------------------------------------------------------------------------------------------------+
| impropers | impropers     | GSD, HOOMD XML, LAMMPS data, PRMTOP, PSF, TPR                                                         |
+-----------+---------------+-------------------------------------------------------------------------------------------------------+


---------------------------------
Adding or modifying TopologyAttrs
---------------------------------

Each of the 


.. _topology-objects:

Topology objects
================

MDAnalysis defines four types of :class:`~MDAnalysis.core.topologyobjects.TopologyObject` by connectivity:

    * :class:`~MDAnalysis.core.topologyobjects.Bond`
    * :class:`~MDAnalysis.core.topologyobjects.Angle`
    * :class:`~MDAnalysis.core.topologyobjects.Dihedral`
    * :class:`~MDAnalysis.core.topologyobjects.Improper`

The value of each topology object can be calculated with :func:`value`, or the name of the topology object (:func:`angle` in this case). See :ref:`topology-objects` for more information. 

Each :class:`~MDAnalysis.core.topologyobjects.TopologyObject` also contain the following attributes:

    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.atoms` : the ordered atoms in the object
    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.indices` : the ordered atom indices in the object
    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.type` : this is either the 'type' of the bond/angle/dihedral/improper, or a tuple of the atom types.
    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.is_guessed` : MDAnalysis can guess bonds. This property records if the object was read from a file or guessed.


Groups of these are held in :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`\ s. The master groups of TopologyObjects are :ref:`accessible as properties of a Universe <universe-properties-label>`. TopologyObjects are typically read from a file with connectivity information (:ref:`see the supported formats here <topologyobject-attr-label>`). However, they can be created in two ways: by adding them to a Universe, or by creating them from an AtomGroup. While the first method adds new objects to the respective Universe master :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`, the second does not.

Bonds can be guessed based on distance and Van der Waals' radii with :meth:`AtomGroup.guess_bonds <MDAnalysis.core.groups.AtomGroup.guess_bonds>`.


Users can also define new TopologyObjects through an :class:`~MDAnalysis.core.groups.AtomGroup`.

--------------------
Adding to a Universe
--------------------



--------------------------
Creating with an AtomGroup
--------------------------

An :class:`~MDAnalysis.core.groups.AtomGroup` can be represented as a bond, angle, dihedral angle, or improper angle :class:`~MDAnalysis.core.topologyobjects.TopologyObject` through the respective properties:

    * :attr:`~MDAnalysis.core.groups.AtomGroup.bond`
    * :attr:`~MDAnalysis.core.groups.AtomGroup.angle`
    * :attr:`~MDAnalysis.core.groups.AtomGroup.dihedral`
    * :attr:`~MDAnalysis.core.groups.AtomGroup.improper`

The :class:`~MDAnalysis.core.groups.AtomGroup` must contain the corresponding number of atoms, in the desired order. For example, a bond cannot be created from three atoms.

.. code-block:: python

    >>> u.atoms[[3, 4, 2]].bond
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "MDAnalysis/package/MDAnalysis/core/groups.py", line 2954, in bond
        "bond only makes sense for a group with exactly 2 atoms")
    ValueError: bond only makes sense for a group with exactly 2 atoms

However, the angle Atom 2 ----- Atom 4 ------ Atom 3 can be calculated, even if the atoms are not connected with bonds.

.. code-block:: python

    >>> a = u.atoms[[3, 4, 2]].angle
    >>> a
    <Angle between: Atom 2, Atom 4, Atom 3>
    >>> a.angle
    <bound method Angle.angle of <Angle between: Atom 2, Atom 4, Atom 3>>
    >>> a.angle()
    47.63986538582528
    >>> a.value()
    47.63986538582528
