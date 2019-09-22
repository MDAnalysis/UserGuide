
.. -*- coding: utf-8 -*-

====================
AtomGroup
====================

A Universe contains all particles in the molecular system, which MDAnalysis calls :code:`Atom`\ s. :code:`Atom`\ s are grouped into :code:`AtomGroup`\ s; the 'master' :code:`AtomGroup` of a Universe is accessible at :code:`Universe.atoms`. Other AtomGroup instances are typically created with :func:`Universe.select_atoms()` or by manipulating another :code:`AtomGroup`, e.g. by slicing.

.. code-block:: python

    >>> u.select_atoms('resname ARG')
    <AtomGroup with 312 atoms>

See :ref:`Selections` for more information.

The :code:`AtomGroup` is probably the most important object in MDAnalysis. Virtually everything can be accessed through an :code:`AtomGroup`. 


Creating an AtomGroup
=====================

--------------------
Indexing and slicing
--------------------

An :code:`AtomGroup` can be indexed and sliced like a list:

.. code-block:: python

    >>> print(u.atoms[0])
    <Atom 1: N of type N of resname MET, resid 1 and segid SYSTEM and altLoc >

Slicing returns another :code:`AtomGroup`. The below code returns an :code:`AtomGroup` of every second element from the first to the 6th element, corresponding to indices 0, 2, and 4.

.. code-block:: python

    >>> ag = u.atoms[0:6:2]
    >>> print(ag)
    <AtomGroup [<Atom 1: N of type N of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 3: H2 of type H of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 5: CA of type C of resname MET, resid 1 and segid SYSTEM and altLoc >]>
    >>> ag.indices
    array([0, 2, 4])


MDAnalysis also supports fancy indexing: passing a :code:`numpy.ndarray` or a :code:`list`. 

.. code-block:: python

    >>> indices = [0, 3, -1, 10, 3]
    >>> u.atoms[indices].indices
    array([    0,     3, 47680,    10,     3])


Boolean indexing allows you to pass in an array of :code:`True` or :code:`False` values to create a new :code:`AtomGroup` from another. The array must be the same length as the original :code:`AtomGroup`. This allows you to select atoms on conditions.

.. code-block:: python

    >>> arr = u.atoms.resnames == 'ARG'
    >>> arr
    array([False, False, False, ..., False, False, False])
    >>> u.atoms[arr]
    <AtomGroup with 312 atoms>


----------------------------------------
Group operators, set methods, and split
----------------------------------------

:code:`AtomGroup`\ s can be compared and combined using group operators.

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

Unlike set methods and atom selection language, concatenation and subtraction keep the order of the atoms as well as duplicates.

.. code-block:: python

    >>> ag1 = u.atoms[1:6]
    >>> ag2 = u.atoms[8:3:-1]
    >>> ag = ag1 + ag2
    >>> ag.indices
    array([1, 2, 3, 4, 5, 8, 7, 6, 5, 4])


MDAnalysis supports a number of set methods. Each of these methods treat the groups as sorted sets of unique :code:`Atom`\ s.


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

.. code-block:: python

    >>> ag1 = u.atoms[1:6]
    >>> ag2 = u.atoms[8:3:-1]
    >>> ag = ag1 | ag2
    >>> ag.indices
    array([1, 2, 3, 4, 5, 6, 7, 8])

:func:`AtomGroup.split()` can create :code:`AtomGroup`\ s by splitting another :code:`AtomGroup` by *level*. 

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

:code:`AtomGroup`\ s can be split by atom, residue, molecule, or segment.

-----------------------
Constructing from Atoms
-----------------------

An :code:`AtomGroup` can be created from an iterable of :code:`Atom`\ s:

.. code-block:: python

    >>> atom1 = u.atoms[4]
    >>> atom2 = u.atoms[6]
    >>> atom3 = u.atoms[2]
    >>> ag = mda.AtomGroup([atom1, atom2, atom3])
    >>> print(ag)
    <AtomGroup [<Atom 5: CA of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 7: CB of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 3: H2 of type H of resname MET, resid 1 and segid SYSTEM and altLoc >]>


Or from providing a list of indices and the Universe that the :code:`Atom`\ s belong to:

.. code-block:: python

    >>> ag = mda.AtomGroup([4, 6, 2], u)
    >>> print(ag)
    <AtomGroup [<Atom 5: CA of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 7: CB of type C of resname MET, resid 1 and segid SYSTEM and altLoc >, <Atom 3: H2 of type H of resname MET, resid 1 and segid SYSTEM and altLoc >]>



Methods
=================

----------------
Topology objects
----------------

* angle
* bond
* dihedral
* improper

----------------------
Groupby and accumulate
----------------------

* Groupby
* accumulate

----------------
Analysis methods
----------------

* bbox
* bsphere
* center
* center_of_geometry
* center_of_mass
* centroid
* dimensions
* forces
* guess_bonds
* pack_into_box
* positions
* unwrap
* velocities
* wrap

--------------------
Manipulation methods
--------------------

* rotate
* rotate_by
* transform
* translate




ResidueGroup and SegmentGroup
=============================


For convenience, chemically meaningful groups of Atoms such as a 
Residue or a Segment (typically a whole molecule or all of the solvent) 
also exist as containers, as well as groups of these units (
ResidueGroup, SegmentGroup).

