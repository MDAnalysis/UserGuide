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

The table below lists attributes that are read from supported formats. These can also be :ref:`added to a Universe <add-topologyattrs-label>` created from a file that does not support them. 

.. include:: scripts/topologyattrs.txt

.. _topologyobject-attr-label:

-----------------------------------
Connectivity information
-----------------------------------

MDAnalysis can also read connectivity information, if the file provides it. These become available as :ref:`topology-objects`, which have additional functionality.

.. include:: scripts/connectivityattrs.txt


.. _add-topologyattrs-label:

----------------------------------
Adding TopologyAttrs
----------------------------------

Each of the attributes above can be added to a Universe if it was not available in the file with :meth:`~MDAnalysis.core.universe.Universe.add_TopologyAttr`.

:meth:`~MDAnalysis.core.universe.Universe.add_TopologyAttr` takes two arguments:

    * :code:`topologyattr` : the singular or plural name of a TopologyAttr, *or* a MDAnalysis TopologyAttr object. This must already have been defined as a :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` (see :ref:`custom-topologyattrs-label` for an example of adding a custom topology attribute).
    * :code:`values` (optional) : if :code:`topologyattr` is a string, the values for that attribute. This can be :code:`None` if the attribute has default values defined, e.g. :code:`bfactors`. 

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF
    psf = mda.Universe(PSF)
    hasattr(psf.atoms, 'bfactors')
    psf.add_TopologyAttr('bfactors')
    psf.atoms.bfactors

One way to modify topology attributes is to simply replace them with :meth:`~MDAnalysis.core.universe.Universe.add_TopologyAttr`:

.. ipython:: python

    psf.add_TopologyAttr('bfactors', range(len(psf.atoms)))
    psf.atoms.bfactors

The number of values provided should correspond with the "level" of the attribute. For example, B-factors are atomic-level information. However, residue names and residue ids apply to residues. See a :ref:`table of attribute levels and default values <topologyattr-defaults-label>` for more information.

----------------------------------
Modifying TopologyAttrs
----------------------------------

Existing topology attributes can be directly modified by assigning new values.

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PDB
    pdb = mda.Universe(PDB)
    pdb.residues[:3].resnames
    pdb.residues[:3].resnames = ['RES1', 'RES2', 'RES3']
    pdb.residues[:3].atoms.resnames


.. note:: 

    This method cannot be used with connectivity attributes, i.e. bonds, angles, dihedrals, and impropers. 


Similarly to adding topology attributes with :meth:`~MDAnalysis.core.universe.Universe.add_TopologyAttr`, the "level" of the attribute matters. Residue attributes can only be assigned to attributes at the Residue or ResidueGroup level. The same applies to attributes for Atoms and Segments. For example, if we tried to assign resnames to an AtomGroup:

.. ipython:: python
    :okexcept:
    
    pdb.residues[0].atoms.resnames = ['new_name']


.. _topologyattr-defaults-label:

-----------------------------------
Default values and attribute levels
-----------------------------------

Topology information in MDAnalysis is always associated with a level: one of atom, residue, or segment. For example, `indices` is Atom information, `resindices` is Residue information, and `segindices` is Segment information. Many topology attributes also have default values, so that they can be :ref:`added to a Universe without providing explicit values <add-topologyattrs-label>`, and expected types. The table below lists which attributes have default values, what they are, and the information level.

.. include:: scripts/topology_defaults.txt



.. _topology-objects:

Topology objects
================

MDAnalysis defines four types of :class:`~MDAnalysis.core.topologyobjects.TopologyObject` by connectivity:

    * :class:`~MDAnalysis.core.topologyobjects.Bond`
    * :class:`~MDAnalysis.core.topologyobjects.Angle`
    * :class:`~MDAnalysis.core.topologyobjects.Dihedral`
    * :class:`~MDAnalysis.core.topologyobjects.Improper`

The value of each topology object can be calculated with :func:`value`, or the name of the topology object. The only exception to this is :class:`~MDAnalysis.core.topologyobjects.Bond`, where the other name for :meth:`~MDAnalysis.core.topologyobjects.Bond.value` is :meth:`~MDAnalysis.core.topologyobjects.Bond.length`.

Each :class:`~MDAnalysis.core.topologyobjects.TopologyObject` also contains the following attributes:

    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.atoms` : the ordered atoms in the object
    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.indices` : the ordered atom indices in the object
    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.type` : this is either the 'type' of the bond/angle/dihedral/improper, or a tuple of the atom types.
    * :attr:`~MDAnalysis.core.topologyobjects.TopologyObject.is_guessed` : MDAnalysis can guess bonds. This property records if the object was read from a file or guessed.


Groups of these are held in :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`\ s. The master groups of TopologyObjects are :ref:`accessible as properties of a Universe <universe-properties-label>`. TopologyObjects are typically read from a file with connectivity information (:ref:`see the supported formats here <topologyobject-attr-label>`). However, they can be created in two ways: by adding them to a Universe, or by creating them from an :class:`~MDAnalysis.core.groups.AtomGroup`. Bonds can be guessed based on distance and Van der Waals' radii with :meth:`AtomGroup.guess_bonds <MDAnalysis.core.groups.AtomGroup.guess_bonds>`.

--------------------
Adding to a Universe
--------------------

As of version 0.21.0, TopologyObjects have specific methods for adding to a Universe:

    * :meth:`~MDAnalysis.core.universe.Universe.add_Bonds`,
    * :meth:`~MDAnalysis.core.universe.Universe.add_Angles`,
    * :meth:`~MDAnalysis.core.universe.Universe.add_Dihedrals`,
    * :meth:`~MDAnalysis.core.universe.Universe.add_Impropers`

These accept the following values:

    * a `~MDAnalysis.core.topologyobjects.TopologyGroup`
    * an iterable of atom indices
    * an iterable of `~MDAnalysis.core.topologyobjects.TopologyObject`\ s

.. todo:: 

    Provide example when 0.21.0 is released


Prior to version 0.21.0, objects could be added to a Universe with :meth:`~MDAnalysis.core.universe.Universe.add_TopologyAttr`. 

.. ipython:: python

    hasattr(pdb, 'angles')
    pdb.add_TopologyAttr('angles', [(0, 1, 2), (2, 3, 4)])
    pdb.angles


Both of these methods add the new objects to the assocated master :class:`~MDAnalysis.core.topologyobjects.TopologyGroup` in the Universe.


--------------------------
Creating with an AtomGroup
--------------------------

An :class:`~MDAnalysis.core.groups.AtomGroup` can be represented as a bond, angle, dihedral angle, or improper angle :class:`~MDAnalysis.core.topologyobjects.TopologyObject` through the respective properties:

    * :attr:`~MDAnalysis.core.groups.AtomGroup.bond`
    * :attr:`~MDAnalysis.core.groups.AtomGroup.angle`
    * :attr:`~MDAnalysis.core.groups.AtomGroup.dihedral`
    * :attr:`~MDAnalysis.core.groups.AtomGroup.improper`

The :class:`~MDAnalysis.core.groups.AtomGroup` must contain the corresponding number of atoms, in the desired order. For example, a bond cannot be created from three atoms.

.. ipython:: python
    :okexcept:

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PDB
    pdb = mda.Universe(PDB)
    pdb.atoms[[3, 4, 2]].bond


However, the angle Atom 2 ----- Atom 4 ------ Atom 3 can be calculated, even if the atoms are not connected with bonds.

.. ipython:: python

    a = pdb.atoms[[3, 4, 2]].angle
    print(a.value())


