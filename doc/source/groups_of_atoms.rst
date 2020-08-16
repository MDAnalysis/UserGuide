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

Access to various data views via :class:`~MDAnalysis.core.groups.AtomGroup` object can be pretty powerful, but also
needs to be used with caution to avoid accessing data outside the intended selection, as mentioned above. Therefore,
we present two use cases showing commonly used applications.

---------------------------------------
Use case: Sequence of residues by chain
---------------------------------------
Let's assume that we have a PDB file, which contains both ATOM and HETATM record types, for protein and small molecule,
respectively. In order to select only ATOM record types and get a list of residues by chain, one needs to call:

.. code-block:: python

    u = mda.Universe(pdb, format="PDB")
    for seg in u.segments:
        p_seg = seg.atoms.select_atoms("record_type ATOM")
        r_list = p_seg.residues.resnames
        print(r_list)

This can be used later on to check whether the residue names belong, for example, to the list of standard residues.
However, it would be **incorrect** to use:

.. code-block:: python

    selected_atoms = u.select_atoms("record_type ATOM")  # CORRECT SYNTAX, BUT INCORRECT RETURN VALUE!
    residues_wrong = selected_atoms.segments.resnames  # CORRECT SYNTAX, BUT INCORRECT RETURN VALUE!
    print(residues_wrong)  # CORRECT SYNTAX, BUT INCORRECT RETURN VALUE!

because, although it nicely returns list of residues in every chain (list of lists, so no need to additionally
loop through segments), it would return all residues, both for ATOM and HETATM record types, because this construct
means "give me all residue names from segments in which there is at least one of the selected atoms".

----------------------------------------
Use case: Atoms list grouped by residues
----------------------------------------
In another use case, let's assume that we would like to check if all the heavy protein backbone and sidechain atoms
are present in every residue. We can easily prepare the expected atom names list in every standard residue, but
we need to obtain the list of atoms in every residue to compare with (using, for example, a symmetrical difference).
In order to do that, we want to avoid HETATMs, OXT and hydrogen (H*) atoms, and we need to loop over segments (chains)
and residues:

.. code-block:: python

    u = mda.Universe(pdb, format="PDB")
    for seg in u.segments:
        p_seg = seg.atoms.select_atoms("record_type ATOM and not (name H* or name OXT)")
        for p_res in p_seg.residues:
            atoms_in_residue = p_seg.select_atoms(f"resid {p_res.resid} and resname {p_res.resname}").names.tolist()
            print(atoms_in_residue)

However, it would be **incorrect** to use in the body of the first loop:

.. code-block:: python

    atoms_in_residue_wrong = p_seg.residues.names  # CORRECT SYNTAX, BUT INCORRECT RETURN VALUE!
    print(atoms_in_residue_wrong)  # CORRECT SYNTAX, BUT INCORRECT RETURN VALUE!

because, although it nicely returns list of atoms in every residue, it would return all atoms (including hydrogen ones),
because this construct means "give me all atom names from residues in which there is at least one of the selected
atoms".

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

    u.atoms[0].fragment  # first atom of lysozyme

and see which fragments are associated with atoms in a smaller :class:`~MDAnalysis.core.groups.AtomGroup`:

.. ipython:: python

    u.atoms[1959:1961].fragments

.. note::

    :attr:`AtomGroup.fragments <MDAnalysis.core.groups.AtomGroup.fragments>` returns a tuple of fragments with at least one :class:`~MDAnalysis.core.groups.Atom` in the :class:`~MDAnalysis.core.groups.AtomGroup`, not a tuple of fragments where *all* Atoms are in the :class:`~MDAnalysis.core.groups.AtomGroup`.