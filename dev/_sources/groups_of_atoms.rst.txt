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

.. ipython:: python
    :okwarning:

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import TPR, XTC
    u = mda.Universe(TPR, XTC)
    ag = u.atoms.select_atoms('resname ARG and name CA')
    ag

Each of these container groups can be accessed through another. The behaviour of this differs by level. For example, the residues of the ``ag`` are the residues that the atoms of ``ag`` belong to. 

.. ipython:: python

    ag.residues

Accessing the atoms of *those* residues, however, returns *all* the atoms in the residues, not just those originally in ``ag``.

.. ipython:: python

    ag.residues.atoms

The same applies to segments.

.. ipython:: python

    ag[:3].segments.atoms

Similarly, an :class:`~MDAnalysis.core.groups.Atom` has direct knowledge of the :class:`~MDAnalysis.core.groups.Residue` and :class:`~MDAnalysis.core.groups.Segment` it belongs to. Note that an :class:`~MDAnalysis.core.groups.Atom` belongs to *one* :class:`~MDAnalysis.core.groups.Residue` and the residue belongs to *one* :class:`~MDAnalysis.core.groups.Segment`, but a :class:`~MDAnalysis.core.groups.Segment` has multiple residues.

.. code-block:: ipython

    In [9]: a = u.atoms[0]

    In [10]: a.residue
    Out[10]: <Residue LYSH, 0>

    In [11]: a.residue.segment
    Out[11]: <Segment seg_0_Protein_A>

    In [12]: a.residue.segment.residues
    Out[12]: <ResidueGroup with 129 residues>

For information on adding custom Residues or Segments, have a look at :ref:`adding-residue`.

Access to other classes via the :class:`~MDAnalysis.core.groups.AtomGroup` object can be pretty powerful, but also
needs to be used with caution to avoid accessing data outside the intended selection. Therefore, we present two use
cases showing commonly used applications, for which we define :class:`~MDAnalysis.core.universe.Universe`
on a simple extract from the PDB file:

.. ipython:: python
    :okwarning:

    import MDAnalysis as mda
    import io
    pdb = io.StringIO("""
    ATOM    414  N   GLY A 402     -51.919   9.578 -14.287  1.00 68.46           N
    ATOM    415  CA  GLY A 402     -52.405  10.954 -14.168  1.00 68.41           C
    ATOM    416  C   GLY A 402     -51.821  11.946 -15.164  1.00 69.71           C
    ATOM    417  O   GLY A 402     -51.958  13.159 -14.968  1.00 69.61           O
    ATOM    418  H   GLY A 402     -52.551   8.935 -14.743  1.00  0.00           H
    ATOM    419  HA3 GLY A 402     -52.225  11.313 -13.155  1.00  0.00           H
    ATOM    420  HA2 GLY A 402     -53.492  10.960 -14.249  1.00  0.00           H
    TER
    HETATM 1929  N1  XYZ A 900     -40.275  19.399 -28.239  1.00  0.00           N1+
    TER
    ATOM   1029  N   ALA B 122     -25.408  19.612 -13.814  1.00 37.52           N
    ATOM   1030  CA  ALA B 122     -26.529  20.537 -14.038  1.00 37.70           C
    ATOM   1031  C   ALA B 122     -26.386  21.914 -13.374  1.00 45.35           C
    ATOM   1032  O   ALA B 122     -26.885  22.904 -13.918  1.00 48.34           O
    ATOM   1033  CB  ALA B 122     -27.835  19.889 -13.613  1.00 37.94           C
    ATOM   1034  H   ALA B 122     -25.636  18.727 -13.385  1.00  0.00           H
    ATOM   1035  HA  ALA B 122     -26.592  20.707 -15.113  1.00  0.00           H
    ATOM   1036  HB1 ALA B 122     -28.658  20.583 -13.783  1.00  0.00           H
    ATOM   1037  HB2 ALA B 122     -27.998  18.983 -14.196  1.00  0.00           H
    ATOM   1038  HB3 ALA B 122     -27.788  19.635 -12.554  1.00  0.00           H
    ATOM   1039  N   GLY B 123     -25.713  21.969 -12.223  1.00 41.18           N
    ATOM   1040  CA  GLY B 123     -25.550  23.204 -11.460  1.00 41.40           C
    ATOM   1041  C   GLY B 123     -24.309  24.018 -11.745  1.00 45.74           C
    ATOM   1042  O   GLY B 123     -24.349  25.234 -11.601  1.00 46.81           O
    ATOM   1043  H   GLY B 123     -25.290  21.133 -11.845  1.00  0.00           H
    ATOM   1044  HA3 GLY B 123     -25.593  22.976 -10.395  1.00  0.00           H
    ATOM   1045  HA2 GLY B 123     -26.430  23.831 -11.600  1.00  0.00           H
    TER
    """)
    u = mda.Universe(pdb, format="PDB")

Use case: Sequence of residues by segment
-----------------------------------------
In order to select only ATOM record types and get a list of residues by segment, one needs to call:

.. ipython:: python

    residues_by_seg = list()
    for seg in u.segments:
        p_seg = seg.atoms.select_atoms("record_type ATOM")
        residues_by_seg.append(p_seg.residues)

Residue names can be extracted using Python's list comprehensions. As required, HETATM record type lines are not
considered:

.. ipython:: python

    [rg.resnames for rg in residues_by_seg]

Note that accessing residues by first selecting the segments of an :class:`~MDAnalysis.core.groups.AtomGroup` returns
all the residues in that segment for both the ATOM and HETATM record types (no memory of the original selection).
The meaning of this is: "give me all residue names from segments in which there is at least one of the selected atoms".

.. ipython:: python

    selected_atoms = u.select_atoms("record_type ATOM")
    all_residues = selected_atoms.segments.residues
    all_residues.resnames

Use case: Atoms list grouped by residues
----------------------------------------
In order to list all the heavy protein backbone and sidechain atoms in every residue, one needs to call:

.. ipython:: python
    :okwarning:

    atoms_in_residues = list()
    for seg in u.segments:
        p_seg = seg.atoms.select_atoms("record_type ATOM and not name H*")
        for p_res in p_seg.residues:
            atoms_in_residues.append(p_seg.select_atoms(f"resid {p_res.resid} and resname {p_res.resname}"))

Atom names can be extracted using Python's list comprehensions. As required, HETATM record type lines and hydrogen
atoms are not considered:

.. ipython:: python

    [ag.names for ag in atoms_in_residues]

The Python syntax can be further simplified by using ``split()`` function:

.. ipython:: python

    rds = u.select_atoms("record_type ATOM and not name H*").split("residue")
    [ag.names for ag in rds]

Note that accessing atoms by first selecting the residues of an :class:`~MDAnalysis.core.groups.AtomGroup` also
returns hydrogen atoms (no memory of the original selection). The meaning of this is "give me all atom names from
residues in which there is at least one of the selected atoms". However, it doesn't contain a nitrogen atom from XYZ
residue as no atoms from this residue were in the :class:`~MDAnalysis.core.groups.AtomGroup`.

.. ipython:: python

    all_atoms_in_residues = list()
    for seg in u.segments:
        p_seg = seg.atoms.select_atoms("record_type ATOM and not name H*")
        all_atoms_in_residues.append(p_seg.residues)
    [atom.names for atom in all_atoms_in_residues]

---------------------------
Fragments
---------------------------

Certain analysis methods in MDAnalysis also make use of additional ways to group atoms. A key concept is a fragment. A fragment is what is typically considered a molecule: an AtomGroup where any atom is reachable from any other atom in the AtomGroup by traversing bonds, and none of its atoms is bonded to any atoms outside the AtomGroup. (A 'molecule' in MDAnalysis methods :ref:`refers to a GROMACS-specific concept <molecule>`). The fragments of a Universe are determined by MDAnalysis as a derived quantity. They can only be determined if bond information is available.

The fragments of an :class:`~MDAnalysis.core.groups.AtomGroup` are accessible via the :attr:`fragments` property. Below is a Universe from a GROMACS TPR file of lysozyme (`PDB ID: 2LYZ <http://www.rcsb.org/structure/2LYZ>`_) with 101 water molecules. While it has 230 residues, there are only 102 fragments: 1 protein and 101 water fragments.

.. code-block:: ipython

    In [12]: len(u.residues)
    Out[12]: 230

    In [13]: len(u.atoms.fragments)
    Out[13]: 102


See :ref:`topology-objects` for more on bonds and which file formats give MDAnalysis bond information.

You can also look at which fragment a particular :class:`~MDAnalysis.core.groups.Atom` belongs to:

.. ipython:: python
    :okwarning:

    u.atoms[0].fragment  # first atom of lysozyme

and see which fragments are associated with atoms in a smaller :class:`~MDAnalysis.core.groups.AtomGroup`:

.. ipython:: python

    u.atoms[1959:1961].fragments

.. note::

    :attr:`AtomGroup.fragments <MDAnalysis.core.groups.AtomGroup.fragments>` returns a tuple of fragments with at least one :class:`~MDAnalysis.core.groups.Atom` in the :class:`~MDAnalysis.core.groups.AtomGroup`, not a tuple of fragments where *all* Atoms are in the :class:`~MDAnalysis.core.groups.AtomGroup`.
