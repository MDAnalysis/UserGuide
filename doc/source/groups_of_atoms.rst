.. -*- coding: utf-8 -*-
.. _groups-of-atoms:

Groups of atoms
===============

MDAnalysis has a hierarchy of :class:`~MDAnalysis.core.groups.Atom` containers that are used throughout the code.

.. image:: images/classes.png

First and foremost is the :class:`~MDAnalysis.core.groups.AtomGroup`. An :class:`~MDAnalysis.core.groups.AtomGroup` is the primary :class:`~MDAnalysis.core.groups.Atom` container; virtually everything can be accessed through it, as detailed in :ref:`atomgroup`. This includes chemically meaningful groups of :class:`~MDAnalysis.core.groups.Atom`\ s such as a :class:`~MDAnalysis.core.groups.Residue` or a :class:`~MDAnalysis.core.groups.Segment`. 

.. _residues-and-segments:

---------------------
Residues and Segments
---------------------

A :class:`~MDAnalysis.core.groups.Residue` is composed of :class:`~MDAnalysis.core.groups.Atom`\ s, and a :class:`~MDAnalysis.core.groups.Segment` is composed of :class:`Residues <MDAnalysis.core.groups.Residue>`.

The corresponding container groups are :class:`~MDAnalysis.core.groups.ResidueGroup` and :class:`~MDAnalysis.core.groups.SegmentGroup`. These have similar properties and available methods as :class:`~MDAnalysis.core.groups.AtomGroup`.

Each of these container groups can be accessed through another. For example:

.. code-block::

    >>> ag = u.atoms.select_atoms('resname ARG and name CA')
    >>> ag
    <AtomGroup with 13 atoms>
    >>> ag.residues
    <ResidueGroup with 13 residues>
    >>> ag.residues.atoms
    <AtomGroup with 312 atoms>
    >>> ag.segments
    <SegmentGroup with 1 segment>
    >>> ag.segments.atoms
    <AtomGroup with 3341 atoms>

Similarly, an :class:`~MDAnalysis.core.groups.Atom` has direct knowledge of the :class:`~MDAnalysis.core.groups.Residue` and :class:`~MDAnalysis.core.groups.Segment` it belongs to.

.. code-block::

    >>> a = u.atoms[0]
    >>> a.residue
    <Residue MET, 1>
    >>> a.residue.segment
    <Segment 4AKE>
    >>> a.residue.segment.residues
    <ResidueGroup with 214 residues>

For information on adding custom Residues or Segments, have a look at :ref:`adding-residue-label`.

---------------------------
Fragments
---------------------------

Certain analysis methods in MDAnalysis also make use of additional ways to group atoms. A key concept is a fragment. A fragment is what is typically considered a molecule: a group of atoms where each atom is bonded to at least one other atom in the fragment, and are not bonded to any atoms outside the fragment. (A 'molecule' in MDAnalysis methods :ref:`refers to a GROMACS-specific concept <molecule-label>`). The fragments of a Universe are determined by MDAnalysis as a derived quantity. They can only be determined if bond information is available.

The fragments of an :class:`~MDAnalysis.core.groups.AtomGroup` are accessible via the :attr:`fragments` property. Below is a Universe from a GROMACS TPR file of lysozyme (`PDB ID: 2LYZ <http://www.rcsb.org/structure/2LYZ>`_) with 101 water molecules. While it has 2263 residues, there are only 102 fragments: 1 protein and 101 water fragments. ::

    >>> import MDAnalysis as mda
    >>> from MDAnalysis.tests.datafiles import TPR2019B3
    >>> u = mda.Universe(TPR209B3)
    >>> len(u.residues)
    2263
    >>> u.residues.resnames
    array(['LYSH', 'VAL', 'PHE', ..., 'SOL', 'SOL', 'SOL'], dtype=object)
    >>> u.atoms.fragments
    (<AtomGroup with 1960 atoms>, <AtomGroup with 3 atoms>, <AtomGroup with 3 atoms>, ..., <AtomGroup with 3 atoms>, <AtomGroup with 3 atoms>, <AtomGroup with 3 atoms>)
    >>> len(u.atoms.fragments)
    102

See :ref:`topology-objects` for more on bonds and which file formats give MDAnalysis bond information.

You can also look at which fragment a particular :class:`~MDAnalysis.core.groups.Atom` belongs to::

    >>> u.atoms[0].fragment  # first atom of lysozyme
    <AtomGroup with 1960 atoms>

and see which fragments are associated with atoms in a smaller :class:`~MDAnalysis.core.groups.AtomGroup`::

    >>> u.atoms[1959:1961].fragments
    (<AtomGroup with 1960 atoms>, <AtomGroup with 3 atoms>)

.. note::

    :attr:`AtomGroup.fragments <MDAnalysis.core.groups.AtomGroup.fragments>` returns a tuple of fragments with at least one :class:`~MDAnalysis.core.groups.Atom` in the :class:`~MDAnalysis.core.groups.AtomGroup`, not a tuple of fragments where *all* Atoms are in the :class:`~MDAnalysis.core.groups.AtomGroup`.