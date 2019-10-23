.. -*- coding: utf-8 -*-
.. _selections:

=======================
Atom selection language
=======================

AtomGroups can be created by selecting atoms using the MDAnalysis atom selection language:

  >>> kalp = u.select_atoms("segid KALP")


The :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` method of a
:class:`~MDAnalysis.core.groups.AtomGroup` or a
:class:`~MDAnalysis.core.universe.Universe` returns an
:class:`~MDAnalysis.core.groups.AtomGroup`. Selections always return an
:class:`~MDAnalysis.core.groups.AtomGroup` with atoms sorted according to their
index in the topology. This is to ensure that there are not any duplicates,
which can happen with complicated selections.

Group subselections using parentheses::

    >>> u.select_atoms("segid DMPC and not (name H* or type OW)")
    <AtomGroup with 3420 atoms>

Pass in :ref:`named AtomGroups as arguments <preexisting-selections>`::

    >>> sph_6 = u.select_atoms("sphzone 6 protein")
    >>> u.select_atoms("around 3 group sph_6", sph_6=sph_6)
    <AtomGroup with 120 atoms>

Turn off :ref:`periodic boundary conditions for geometric keywords <geometric>` with ``periodic=False``::

    >>> u.select_atoms("around 6 protein", periodic=False)

Create :ref:`dynamic UpdatingAtomGroups <dynamic-selections>` with ``updating=True``::

    >>> u.select_atoms("prop x < 5 and prop y < 5 and prop z < 5", updating=True)
    <AtomGroup with 9 atoms, with selection 'prop x < 5 and prop y < 5 and prop z < 5' on the entire Universe.>

It is possible to export selections for external software
packages with the help of :ref:`selection-exporters`.

Selection Keywords
==================

The following describes all selection keywords currently understood by the
selection parser. The following applies to all selections:

* Keywords are case sensitive.
* Atoms are automatically sequentially ordered in a resulting selection (see
  notes below on :ref:`ordered-selections-label` for how to circumvent this if
  necessary).
* Selections are parsed left to right and parentheses can be used for
  grouping.
* Currently, only "stemming" is implemented as a primitive form of pattern
  matching: Using the ``*`` character in a string such as ``GL*`` selects
  all strings that start with "GL" such as "GLU", "GLY", "GLX29", "GLN". Only terminal wildcards (i.e. matching the last part of a name) are currently supported. 

.. note::

    MDAnalysis will ignore everything after the ``*``. ``u.select_atoms("resname *E")`` will not select atoms whose residue name ends in E, but instead select every atom.


Simple selections
-----------------

protein
    Selects atoms that belong to a :ref:`hard-coded set of standard protein residue names <protein-selection>`.

backbone
    Selects the backbone atoms of a hard-coded set of protein residues. These atoms have the names: CA, C, O, N.

nucleic
    Selects atoms that belong to a :ref:`hard-coded set of standard nucleic residue names <nucleic-selection>`.

nucleicbackbone
    Selects the backbone atoms of a hard-coded set of nucleic residues. These atoms have the names: P, O5', C5', C3', O3'

nucleicbase
    Selects the atoms in :ref:`nucleobases <nucleobase-selection>`.

nucleicsugar
    Selects the atoms in nucleic sugars. These have the names: C1', C2', C3', C4', O2', O4', O3'

segid *seg-name*
    select by segid (as given in the topology), e.g. ``segid 4AKE`` or
    ``segid DMPC``

resid *residue-number-range*
    resid can take a single residue number or a range of numbers. A range
    consists of two numbers separated by a colon (inclusive) such
    as ``resid 1:5``. A residue number ("resid") is taken directly from the
    topology.

resnum *residue-number-range*
    resnum is an alias of resid.

resname *residue-name*
    select by residue name, e.g. ``resname LYS``

name *atom-name*
    select by atom name (as given in the topology). Often, this is force
    field dependent. Example: ``name CA`` (for C-alpha atoms) or ``name
    OW`` (for SPC water oxygen)

type *atom-type*
    select by atom type; this is either a string or a number and depends on
    the force field; it is read from the topology file (e.g. the CHARMM PSF
    file contains numeric atom types). This uses the ``Atom.type`` :ref:`topology attribute <topology-attributes>`.

atom *seg-name residue-number atom-name*
    a selector for a single atom consisting of segid resid atomname,
    e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of the
    DMPC segment

altloc *alternative-location*
    a selection for atoms where alternative locations are available, which is
    often the case with high-resolution crystal structures
    e.g. :code:`resid 4 and resname ALA and altloc B` selects only the atoms of ALA-4
    that have an altloc B record.

moltype *molecule-type*
    select by the ``moltype`` :ref:`topology attribute <topology-attributes>`, e.g. ``moltype Protein_A``. At the moment, only the TPR format defines the ``moltype``.

Boolean
-------

not
    all atoms not in the selection, e.g. ``not protein`` selects all atoms
    that aren't part of a protein

and
    the intersection of two selections, i.e. the boolean and. e.g. ``protein and not resname ALA`` selects all atoms that belong to a protein but are not in an alanine residue

or
    the union of two selections, i.e. the boolean or. e.g. ``protein and not (resname ALA or resname LYS)`` selects all atoms that belong to a protein, but are not in a lysine or alanine residue

.. _geometric:

Geometric
---------

The geometric keywords below all implement periodic boundary conditions by default when valid cell dimensions are accessible from the Universe. This can be turned off by passing in the keyword ``periodic=False``::

    u.select_atoms("around 6 protein", periodic=False)

around *distance selection*
    selects all atoms a certain cutoff away from another selection,
    e.g. ``around 3.5 protein`` selects all atoms not belonging to protein
    that are within 3.5 Angstroms from the protein

sphzone *externalRadius selection*
    selects all atoms within a spherical zone centered in the center of
    geometry (COG) of a given selection, e.g. ``sphzone 6.0 ( protein and (
    resid 130 or resid 80 ) )`` selects the center of geometry of protein,
    resid 130, resid 80 and creates a sphere of radius 6.0 around the COG.

sphlayer *innerRadius externalRadius selection*
    selects all atoms within a spherical layer centered in the center of
    geometry (COG) of a given selection, e.g., ``sphlayer 2.4 6.0 ( protein
    and ( resid 130 or resid 80 ) )`` selects the center of geometry of
    protein, resid 130, resid 80 and creates a spherical layer of inner
    radius 2.4 and external radius 6.0 around the COG.

cyzone *externalRadius zMax zMin selection*
    selects all atoms within a cylindric zone centered in the center of
    geometry (COG) of a given selection, e.g. ``cyzone 15 4 -8 protein and
    resid 42`` selects the center of geometry of protein and resid 42, and
    creates a cylinder of external radius 15 centered on the COG. In z, the
    cylinder extends from 4 above the COG to 8 below. Positive values for
    *zMin*, or negative ones for *zMax*, are allowed.

cylayer *innerRadius externalRadius zMax zMin selection*
    selects all atoms within a cylindric layer centered in the center of
    geometry (COG) of a given selection, e.g. ``cylayer 5 10 10 -8
    protein`` selects the center of geometry of protein, and creates a
    cylindrical layer of inner radius 5, external radius 10 centered on the
    COG. In z, the cylinder extends from 10 above the COG to 8
    below. Positive values for *zMin*, or negative ones for *zMax*, are
    allowed.

point *x y z distance*
    selects all atoms within a cutoff of a point in space, make sure
    coordinate is separated by spaces, e.g. ``point 5.0 5.0 5.0 3.5``
    selects all atoms within 3.5 Angstroms of the coordinate (5.0, 5.0,
    5.0)

prop *[abs] property operator value*
    selects atoms based on position, using *property* **x**, **y**, or
    **z** coordinate. Supports the **abs** keyword (for absolute value) and
    the following *operators*: **<, >, <=, >=, ==, !=**. For example,
    ``prop z >= 5.0`` selects all atoms with z coordinate greater than 5.0;
    ``prop abs z <= 5.0`` selects all atoms within -5.0 <= z <= 5.0.


Similarity and connectivity
---------------------------

same *subkeyword* as *selection*
    selects all atoms that have the same *subkeyword* value as any atom in
    *selection*. Allowed *subkeyword* values are the atom properties: ``name,
    type, resname, resid, resnum, segid, mass, charge, radius, bfactor``, the
    groups an atom belong to: ``residue, segment, fragment``, and the atom
    coordinates ``x, y, z``. (Note that ``bfactor`` currently only works for MMTF formats.) e.g. ``same charge as protein`` selects all atoms that have the same charge as any atom in protein.

byres *selection*
    selects all atoms that are in the same segment and residue as selection,
    e.g. specify the subselection after the byres keyword.  ``byres`` is a
    shortcut to ``same residue as``

bonded *selection*
    selects all atoms that are bonded to selection
    e.g.: ``name H and bonded name N`` selects only hydrogens bonded to
    nitrogens

Index
-----
index *index-range*
    selects all atoms within a range of (0-based) inclusive indices,
    e.g. ``index 0`` selects the first atom in the universe; ``index 5:10``
    selects the 6th through 11th atoms, inclusive. This uses the ``Atom.index`` :ref:`topology attribute <topology-attributes>`.

bynum *number-range*
    selects all atoms within a range of (1-based) inclusive indices,
    e.g. ``bynum 1`` selects the first atom in the universe; ``bynum 5:10``
    selects 5th through 10th atoms, inclusive.

    .. note::

        These are **not** the same as the 1-indexed ``Atom.id`` :ref:`topology attribute <topology-attributes>`. ``bynum`` simply adds 1 to the 0-indexed ``Atom.index``.


.. _preexisting-selections:

Preexisting selections and modifiers
------------------------------------

group `group-name`
    selects the atoms in the :class:`AtomGroup` passed to the function as an
    argument named `group-name`. Only the atoms common to `group-name` and the
    instance :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` was called
    from will be considered, unless ``group`` is preceded by the ``global``
    keyword. `group-name` will be included in the parsing just by comparison of
    atom indices. This means that it is up to the user to make sure the
    `group-name` group was defined in an appropriate :class:`Universe`.

global *selection*
    by default, when issuing
    :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from an
    :class:`~MDAnalysis.core.groups.AtomGroup`, selections and subselections
    are returned intersected with the atoms of that instance.  Prefixing a
    selection term with ``global`` causes its selection to be returned in its
    entirety.  As an example, the ``global`` keyword allows for
    ``lipids.select_atoms("around 10 global protein")`` --- where ``lipids`` is
    a group that does not contain any proteins. Were ``global`` absent, the
    result would be an empty selection since the ``protein`` subselection would
    itself be empty.  When calling
    :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from a
    :class:`~MDAnalysis.core.universe.Universe`, ``global`` is ignored.


.. _dynamic-selections:

Dynamic selections
==================

By default :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` returns an
:class:`~MDAnalysis.core.groups.AtomGroup`, in which the list of atoms is
constant across trajectory frame changes. If
:meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` is invoked with named
argument ``updating`` set to ``True``, an
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup` instance will be returned
instead. It behaves just like an :class:`~MDAnalysis.core.groups.AtomGroup`
object, with the difference that the selection expressions are re-evaluated
every time the trajectory frame changes (this happens lazily, only when the
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup` object is accessed so that
there is no redundant updating going on)::

    # A dynamic selection of corner atoms:
    >>> ag_updating = universe.select_atoms("prop x < 5 and prop y < 5 and prop z < 5", updating=True)
    >>> ag_updating
    <UpdatingAtomGroup with 9 atoms>
    >>> universe.trajectory.next()
    >>> ag_updating
    <UpdatingAtomGroup with 14 atoms>

Using the ``group`` selection keyword for
:ref:`preexisting selections < preexisting-selections>`, one can
make updating selections depend on
:class:`~MDAnalysis.core.groups.AtomGroup`, or even other
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup`, instances.
Likewise, making an updating selection from an already updating group will
cause later updates to also reflect the updating of the base group::

    >>> chained_ag_updating = ag_updating.select_atoms("resid 1:1000", updating=True)
    >>> chained_ag_updating
    <UpdatingAtomGroup with 3 atoms>
    >>> universe.trajectory.next()
    >>> chained_ag_updating
    <UpdatingAtomGroup with 7 atoms>

Finally, a non-updating selection or a slicing/addition operation made on an
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup` will return a static
:class:`~MDAnalysis.core.groups.AtomGroup`, which will no longer update
across frames::

    >>> static_ag = ag_updating.select_atoms("resid 1:1000")
    >>> static_ag
    <UpdatingAtomGroup with 3 atoms>
    >>> universe.trajectory.next()
    >>> static_ag
    <UpdatingAtomGroup with 3 atoms>

.. _ordered-selections-label:

Ordered selections
==================

:meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` sorts the atoms
in the :class:`~MDAnalysis.core.groups.AtomGroup` by atom index before
returning them (this is to eliminate possible duplicates in the
selection). If the ordering of atoms is crucial (for instance when
describing angles or dihedrals) or if duplicate atoms are required
then one has to concatenate multiple AtomGroups, which does not sort
them.

The most straightforward way to concatenate two AtomGroups is by using the
``+`` operator::
    >>> ordered = u.select_atoms("resid 3 and name CA") + u.select_atoms("resid 2 and name CA")
    >>> list(ordered)
    [<Atom 46: CA of type C of resname ILE, resid 3 and segid SYSTEM and altLoc >, <Atom 22: CA of type C of resname ARG, resid 2 and segid SYSTEM and altLoc >]

A shortcut is to provide *two or more* selections to
:meth:`~MDAnalysis.core.universe.Universe.select_atoms`, which then
does the concatenation automatically::
    >>> list(u.select_atoms("resid 3 and name CA", "resid 2 and name CA"))
    [<Atom 46: CA of type C of resname ILE, resid 3 and segid SYSTEM and altLoc >, <Atom 22: CA of type C of resname ARG, resid 2 and segid SYSTEM and altLoc >]

Just for comparison to show that a single selection string does not
work as one might expect::
    >>> list(u.select_atoms("(resid 3 or resid 2) and name CA"))
    [<Atom 22: CA of type C of resname ARG, resid 2 and segid SYSTEM and altLoc >, <Atom 46: CA of type C of resname ILE, resid 3 and segid SYSTEM and altLoc >]

