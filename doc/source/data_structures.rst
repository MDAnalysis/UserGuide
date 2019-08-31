.. -*- coding: utf-8 -*-

====================
Data structures
====================


Universes
====================

| *If you wish to make an apple pie from scratch, you must first invent the universe.*
| - Carl Sagan, Cosmos


The Universe class ties a topology and a trajectory together. 
Almost all code in MDAnalysis starts with a Universe.

Words words.

-------------------
Creating a Universe
-------------------

Loading from files
------------------

.. code-block:: python

    u = Universe(topology, trajectory)          # read system from file(s)
    u = Universe(pdbfile)                       # read atoms and coordinates from PDB or GRO
    u = Universe(topology, [traj1, traj2, ...]) # read from a list of trajectories
    u = Universe(topology, traj1, traj2, ...)   # read from multiple trajectories


Many words here

Constructing from AtomGroups
----------------------------

Merge()


Constructing from scratch
-------------------------

Universe.empty()



AtomGroup
====================
AtomGroup instances can be easily created 
(e.g., from an AtomGroup.select_atoms() selection or simply by slicing).



ResidueGroup and SegmentGroup
=============================


For convenience, chemically meaningful groups of Atoms such as a 
Residue or a Segment (typically a whole molecule or all of the solvent) 
also exist as containers, as well as groups of these units (
ResidueGroup, SegmentGroup).




Topology
====================

MDA Topology =/= 'force field' topology. Connectivity only. See ParmEd etc for other stuff.