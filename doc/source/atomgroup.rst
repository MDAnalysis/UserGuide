.. -*- coding: utf-8 -*-
.. _atomgroup::

AtomGroup
====================

A :class:`~MDAnalysis.core.universe.Universe` contains all particles in the molecular system. MDAnalysis calls a particle an :class:`~MDAnalysis.core.groups.Atom`, regardless of whether it really is (e.g. it may be a united-atom particle or coarse-grained bead). :class:`~MDAnalysis.core.groups.Atom`\ s are grouped with an :class:`~MDAnalysis.core.groups.AtomGroup`; the 'master' :class:`~MDAnalysis.core.groups.AtomGroup` of a Universe is accessible at :code:`Universe.atoms`. 

.. note::
    The :class:`~MDAnalysis.core.groups.AtomGroup` is probably the most important object in MDAnalysis. Virtually everything can be accessed through an :class:`~MDAnalysis.core.groups.AtomGroup`. 

-----------------------
Creating an AtomGroup
-----------------------


Atom selection language
-----------------------

AtomGroup instances are typically created with :meth:`Universe.select_atoms <MDAnalysis.core.universe.Universe.select_atoms>` or by manipulating another :class:`~MDAnalysis.core.groups.AtomGroup`, e.g. by slicing.

.. code-block:: python

    >>> u.select_atoms('resname ARG')
    <AtomGroup with 312 atoms>

See :ref:`Selections` for more information.


Indexing and slicing
--------------------

An :class:`~MDAnalysis.core.groups.AtomGroup` can be indexed and sliced like a list:

.. code-block:: python

    >>> print(u.atoms[0])
    <Atom 1: N of type N of resname MET, resid 1 and segid SYSTEM and altLoc >

Slicing returns another :class:`~MDAnalysis.core.groups.AtomGroup`. The below code returns an :class:`~MDAnalysis.core.groups.AtomGroup` of every second element from the first to the 6th element, corresponding to indices 0, 2, and 4.

.. code-block:: python

    >>> ag = u.atoms[0:6:2]
    >>> print(ag)
    <AtomGroup [<Atom 1: N of type N of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 3: H2 of type H of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 5: CA of type C of resname MET, resid 1 and segid SYSTEM and altLoc >]>
    >>> ag.indices
    array([0, 2, 4])


MDAnalysis also supports fancy indexing: passing a :class:`~numpy.ndarray` or a :class:`~list`. 

.. code-block:: python

    >>> indices = [0, 3, -1, 10, 3]
    >>> u.atoms[indices].indices
    array([    0,     3, 47680,    10,     3])


Boolean indexing allows you to pass in an array of :code:`True` or :code:`False` values to create a new :class:`~MDAnalysis.core.groups.AtomGroup` from another. The array must be the same length as the original :class:`~MDAnalysis.core.groups.AtomGroup`. This allows you to select atoms on conditions.

.. code-block:: python

    >>> arr = u.atoms.resnames == 'ARG'
    >>> arr
    array([False, False, False, ..., False, False, False])
    >>> u.atoms[arr]
    <AtomGroup with 312 atoms>


Group operators and set methods
-------------------------------

MDAnalysis supports a number of ways to compare :class:`~MDAnalysis.core.groups.AtomGroup`\ s or construct a new one: group operators (e.g. :meth:`~MDAnalysis.core.groups.AtomGroup.concatenate`, :meth:`~MDAnalysis.core.groups.AtomGroup.subtract`) and set methods (e.g. :meth:`~MDAnalysis.core.groups.AtomGroup.union`, :meth:`~MDAnalysis.core.groups.AtomGroup.difference`). Group operators achieve a similar outcome to set methods. However, a key difference is that :meth:`~MDAnalysis.core.groups.AtomGroup.concatenate` and :meth:`~MDAnalysis.core.groups.AtomGroup.subtract` preserve the order of the atoms and any duplicates. :meth:`~MDAnalysis.core.groups.AtomGroup.union` and :meth:`~MDAnalysis.core.groups.AtomGroup.difference` return an :class:`~MDAnalysis.core.groups.AtomGroup` where each atom is unique, and ordered by its topology index. 

.. code-block:: python

    >>> ag1 = u.atoms[1:6]
    >>> ag2 = u.atoms[8:3:-1]
    >>> concat = ag1 + ag2
    >>> concat.indices
    array([1, 2, 3, 4, 5, 8, 7, 6, 5, 4])
    >>> union = ag1 | ag2
    >>> union.indices
    array([1, 2, 3, 4, 5, 6, 7, 8])


**Available operators**

Unlike set methods and atom selection language, concatenation and subtraction keep the order of the atoms as well as duplicates.

+-------------------------------+------------+----------------------------+
| Operation                     | Equivalent | Result                     |
+===============================+============+============================+
| ``len(s)``                    |            | number of atoms            |
|                               |            | in the group               |
+-------------------------------+------------+----------------------------+
| ``s == t``                    |            | test if ``s`` and ``t``    |
|                               |            | contain the same elements  |
|                               |            | in the same order          |
+-------------------------------+------------+----------------------------+
| ``s.concatenate(t)``          | ``s + t``  | new Group with elements    |
|                               |            | from ``s`` and from ``t``  |
+-------------------------------+------------+----------------------------+
| ``s.subtract(t)``             |            | new Group with elements    |
|                               |            | from ``s`` that are not    |
|                               |            | in ``t``                   |
+-------------------------------+------------+----------------------------+

**Available set methods**

Each of these methods create groups that are sorted sets of unique :class:`~MDAnalysis.core.groups.Atom`\ s.

+-------------------------------+------------+----------------------------+
| Operation                     | Equivalent | Result                     |
+===============================+============+============================+
| ``s.isdisjoint(t)``           |            | ``True`` if ``s`` and      |
|                               |            | ``t`` do not share         |
|                               |            | elements                   |
+-------------------------------+------------+----------------------------+
| ``s.issubset(t)``             |            | test if all elements of    |
|                               |            | ``s`` are part of ``t``    |
+-------------------------------+------------+----------------------------+
| ``s.is_strict_subset(t)``     |            | test if all elements of    |
|                               |            | ``s`` are part of ``t``,   |
|                               |            | and ``s != t``             |
+-------------------------------+------------+----------------------------+
| ``s.issuperset(t)``           |            | test if all elements of    |
|                               |            | ``t`` are part of ``s``    |
+-------------------------------+------------+----------------------------+
| ``s.is_strict_superset(t)``   |            | test if all elements of    |
|                               |            | ``t`` are part of ``s``,   |
|                               |            | and ``s != t``             |
+-------------------------------+------------+----------------------------+
| ``s.union(t)``                | ``s | t``  | new Group with elements    |
|                               |            | from both ``s`` and ``t``  |
+-------------------------------+------------+----------------------------+
| ``s.intersection(t)``         | ``s & t``  | new Group with elements    |
|                               |            | common to ``s`` and ``t``  |
+-------------------------------+------------+----------------------------+
| ``s.difference(t)``           | ``s - t``  | new Group with elements of |
|                               |            | ``s`` that are not in ``t``|
+-------------------------------+------------+----------------------------+
| ``s.symmetric_difference(t)`` | ``s ^ t``  | new Group with elements    |
|                               |            | that are part of ``s`` or  |
|                               |            | ``t`` but not both         |
+-------------------------------+------------+----------------------------+

Groupby and split
-----------------

An :class:`~MDAnalysis.core.groups.AtomGroup` can be constructed from another by separating atoms by properties. 

:meth:`AtomGroup.split <MDAnalysis.core.groups.AtomGroup.split>` can create a list of :class:`~MDAnalysis.core.groups.AtomGroup`\ s by splitting another :class:`~MDAnalysis.core.groups.AtomGroup` by the 'level' of connectivity: one of *atom*, *residue*, *molecule*, or *segment*. 

.. code-block:: python

    >>> ag1 = u.atoms[:100]
    >>> ag1
    <AtomGroup with 100 atoms>
    >>> ag1.split('residue')
    [<AtomGroup with 19 atoms>,
    <AtomGroup with 24 atoms>,
    <AtomGroup with 19 atoms>,
    <AtomGroup with 19 atoms>,
    <AtomGroup with 19 atoms>]


An :class:`~MDAnalysis.core.groups.AtomGroup` can also be separated according to values of `topology attributes <../topology_system.html#topology-attributes>`_ to produce a dictionary of :code:`{value:AtomGroup}`. 

.. code-block::

    >>> u = mda.Universe(PSF, DCD)
    >>> u.atoms.groupby('masses')
    {32.06: <AtomGroup with 7 atoms>, 1.008: <AtomGroup with 1685 atoms>, 12.011: <AtomGroup with 1040 atoms>, 14.007: <AtomGroup with 289 atoms>, 15.999: <AtomGroup with 320 atoms>}

Passing in multiple attributes groups them in order:

.. code-block::

    >>> u.atoms.groupby(['masses', 'resnames'])
    {(32.06, 'MET'): <AtomGroup with 6 atoms>, (32.06, 'CYS'): <AtomGroup with 1 atom>, (1.008, 'LEU'): <AtomGroup with 176 atoms>, (1.008, 'SER'): <AtomGroup with 25 atoms>, (1.008, 'TYR'): <AtomGroup with 63 atoms>, (1.008, 'ARG'): <AtomGroup with 169 atoms>, (1.008, 'GLU'): <AtomGroup with 108 atoms>, (1.008, 'GLY'): <AtomGroup with 60 atoms>, (1.008, 'ASN'): <AtomGroup with 24 atoms>, ..., }

Constructing from Atoms
-----------------------

An :class:`~MDAnalysis.core.groups.AtomGroup` can be created from an iterable of :class:`~MDAnalysis.core.groups.Atom`\ s:

.. code-block:: python

    >>> atom1 = u.atoms[4]
    >>> atom2 = u.atoms[6]
    >>> atom3 = u.atoms[2]
    >>> ag = mda.AtomGroup([atom1, atom2, atom3])
    >>> print(ag)
    <AtomGroup [<Atom 5: CA of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 7: CB of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 3: H2 of type H of resname MET, resid 1 and segid SYSTEM and altLoc >]>


Or from providing a list of indices and the Universe that the :class:`~MDAnalysis.core.groups.Atom`\ s belong to:

.. code-block:: python

    >>> ag = mda.AtomGroup([4, 6, 2], u)
    >>> print(ag)
    <AtomGroup [<Atom 5: CA of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 7: CB of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 3: H2 of type H of resname MET, resid 1 and segid SYSTEM and altLoc >]>

Order and uniqueness
-----------------------

These methods of creating an :class:`~MDAnalysis.core.groups.AtomGroup` result in a sorted, unique list of atoms:

    * Atom selection language
    * Slicing
    * Boolean indexing
    * Set methods
    * :meth:`AtomGroup.split <MDAnalysis.core.groups.AtomGroup.split>` and :meth:`AtomGroup.groupby <MDAnalysis.core.groups.AtomGroup.groupby>`
    
These methods return a user-ordered :class:`~MDAnalysis.core.groups.AtomGroup` that can contain duplicates:

    * Fancy indexing (with arrays or lists)
    * Group operations (:meth:`AtomGroup.concatenate <MDAnalysis.core.groups.AtomGroup.concatenate>` and :meth:`AtomGroup.subtract <MDAnalysis.core.groups.AtomGroup.subtract>`)
    * Constructing directly from :class:`~MDAnalysis.core.groups.Atom`\ s

Empty AtomGroups
----------------

Empty AtomGroups can be constructed from each method. For example, using selection language::

    >>> u.select_atoms("resname DOES_NOT_EXIST")
    <AtomGroup with 0 atoms>

and indexing::

    >>> u.atoms[6:6]
    <AtomGroup with 0 atoms>

or set methods::

    >>> u.atoms - u.atoms
    <AtomGroup with 0 atoms>

Creating an :class:`~MDAnalysis.core.groups.AtomGroup` from an empty list requires passing in a :class:`~MDAnalysis.core.universe.Universe`. ::

    >>> ag = mda.AtomGroup([], u)
    >>> ag
    <AtomGroup with 0 atoms>

Empty AtomGroups have a length of 0 and evaluate to :code:`False` in a boolean context. ::

    >>> null = u.atoms[[]]
    >>> len(null)
    0
    >>> bool(null)
    False

AtomGroups with Atoms evaluate to :code:`True` in a boolean context::

    >>> ag = u.atoms[:5]
    >>> len(ag)
    5
    >>> bool(ag)
    True

-------
Methods
-------

Most of the analysis functionality in MDAnalysis is implemented in :mod:`the analysis module <analysis>`, but many interesting methods can be accessed from an :class:`~MDAnalysis.core.groups.AtomGroup` directly. 

Topology objects
----------------

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

The value of each topology object can be calculated with :func:`value`, or the name of the topology object (:func:`angle` in this case). See :ref:`topology-objects` for more information. 


Topology-specific methods
-------------------------

A number of analysis and transformation methods are defined for :class:`~MDAnalysis.core.groups.AtomGroup` that require specific properties to be available. The primary requirement is usually coordinates, or the `positions` attribute. With coordinates, you can easily compute a center of geometry:

.. code-block::

    >>> u.atoms.center_of_geometry()
    array([-0.04223882,  0.01418196, -0.03504874])

When atom masses are defined, a number of other methods become available. These include:

* :func:`~MDAnalysis.core.groups.AtomGroup.center_of_mass`
* :func:`~MDAnalysis.core.groups.AtomGroup.radius_of_gyration`
* :func:`~MDAnalysis.core.groups.AtomGroup.moment_of_inertia`
* :func:`~MDAnalysis.core.groups.AtomGroup.asphericity`
* :func:`~MDAnalysis.core.groups.AtomGroup.align_principal_axes`

See :ref:`the Topology system <topology-system>` for more information about which attributes are required for which :class:`~MDAnalysis.core.groups.AtomGroup` methods.
