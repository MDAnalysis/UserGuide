.. -*- coding: utf-8 -*-

====================
Universes
====================

| *If you wish to make an apple pie from scratch, you must first invent the universe.*
| - Carl Sagan, Cosmos

| *In the beginning the Universe was created.*
| *This has made a lot of people very angry and been widely regarded as a bad move.*
| ― Douglas Adams, The Restaurant at the End of the Universe

The Universe class ties a topology and a trajectory together. 
Almost all code in MDAnalysis starts with a Universe.

Creating a Universe
===================

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

