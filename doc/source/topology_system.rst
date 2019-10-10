.. -*- coding: utf-8 -*-

=====================
The topology system
=====================

MDAnalysis groups static data about a :code:`Universe` into its topology. This is typically loaded from a topology file. Topology information falls into 3 categories:

    * Atom containers (Residues and Segments)
    * Atom attributes (e.g. name, mass, bfactor)
    * Topology objects: bonds, angles, dihedrals, impropers

Users will almost never interact directly with a :code:`Topology`. Modifying atom containers or topology attributes is typically done through :code:`Universe`. Methods for viewing containers or topology attributes, or for calculating topology object values, are accessed through :code:`AtomGroup`.

Atom containers
===============

To add a :code:`Residue` or :code:`Segment` to a topology, use the :code:`Universe.add_Residue` or :code:`Universe.add_Segment` methods.

.. code-block::

    >>> u = mda.Universe(PSF, DCD)
    >>> u.segments
    <SegmentGroup with 1 segment>
    >>> u.segments.segids
    array(['4AKE'], dtype=object)
    >>> newseg = u.add_Segment(segid='X')
    >>> u.segments.segids
    array(['4AKE', 'X'], dtype=object)
    >>> newseg.atoms
    <AtomGroup with 0 atoms>

To assign the last 100 residues from the :code:`Universe` to this new :code:`Segment`:

.. code-block::

    >>> u.residues[-100:].segments = newseg
    >>> newseg.atoms
    <AtomGroup with 1600 atoms>


Topology attributes
===================

MDAnalysis supports a range of topology attributes for each :code:`Atom`. These are only available if read from the files loaded into the :code:`Universe`. 

Several topology attributes are central to MDAnalysis, and are derived for each Universe regardless of format:






---------------------------------
Adding or modifying TopologyAttrs
---------------------------------

See this notebook on colouring a protein by RMSF with the tempfactor attribute.

Topology objects
================

MDAnalysis defines four :code:`TopologyObject`\ s by connectivity:

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

Bonds can be guessed based on distance and Van der Waals' radii with :code:`guess_bonds`:


Users can also define new TopologyObjects through an :code:`AtomGroup`.



----------------------------
Calculating values
----------------------------

