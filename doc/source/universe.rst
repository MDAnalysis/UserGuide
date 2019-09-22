.. -*- coding: utf-8 -*-

====================
Universe
====================

    *If you wish to make an apple pie from scratch, you must first invent the universe.*

    -- Carl Sagan, Cosmos


:code:`MDAnalysis.Universe` contains all the information describing a molecular dynamics system. Almost all code in MDAnalysis begins with a Universe. It has two key properties:

* :code:`atoms`: an :code:`AtomGroup` of the system's atoms, providing access to important analysis methods (described below)
* :code:`trajectory`: the currently loaded trajectory reader


:code:`MDAnalysis.Universe` ties the static information from the "topology" (e.g. atom identities) to dynamically updating information from the "trajectory" (e.g. coordinates). A cornerstone of MDAnalysis' ontology is that an entire trajectory is never loaded into memory. Instead, the :code:`trajectory` attribute provides a view on a specific frame of the trajectory. This allows the analysis of arbitrarily long trajectories without a significant impact on memory. 

Creating a Universe
===================

------------------
Loading from files
------------------

A Universe is typically created from a "topology" file, with optional "trajectory" file/s. Trajectory files must have the coordinates in the same order as atoms in the topology. See FORMATS for the topology and trajectory formats supported by MDAnalysis.

.. code-block:: python

    u = Universe(topology, trajectory)          
    u = Universe(pdbfile)                       # read atoms and coordinates from PDB or GRO
    u = Universe(topology, [traj1, traj2, ...]) # read from a list of trajectories
    u = Universe(topology, traj1, traj2, ...)   # read from multiple trajectories

The line between topology and trajectory files is quite blurry. For example, a PDB or GRO file is considered both a topology and a trajectory file. The difference is that a **topology file** provides static information, such as atom identities (name, mass, etc.), charges, and bond connectivity. A **trajectory file** provides dynamic information, such as coordinates, velocities, forces, and box dimensions. 

If only a single file is provided, MDAnalysis tries to read both topology and trajectory information from it. If multiple trajectory files are provided, coordinates are loaded in the order given. If :code:`continuous=True` is passed, the time steps are modified for a time-continuous trajectory.

The default arguments should create a Universe suited for most analysis applications. However, the :code:`Universe` constructor takes a number of optional arguments.


**The following options specify how to treat the input:**

* :code:`format`: the file format of the trajectory file/s. (default: None, formats are guessed)
* :code:`topology_format`: the file format of the topology file. (default: None, formats are guessed)
* :code:`all_coordinates`: whether to read coordinate information from the first file (default: False. Ignored when only one file is provided)
* :code:`continuous`: whether to give multiple trajectory files continuous time steps. (default: False.)

.. code-block:: python

    >>> import MDAnalysis as mda
    >>> from MDAnalysis.tests.datafiles import PDB, GRO, XTC
    >>> u1 = mda.Universe(GRO, XTC, XTC, all_coordinates=True)
    >>> u1.trajectory
    <ChainReader containing adk_oplsaa.xtc, adk_oplsaa.xtc with 21 frames of 47681 atoms>
    >>> u2 = mda.Universe(GRO, XTC, XTC, all_coordinates=False, continuous=False)
    >>> u2.trajectory
    <ChainReader containing adk_oplsaa.xtc, adk_oplsaa.xtc with 20 frames of 47681 atoms>
    >>> print([int(ts.time) for ts in u2.trajectory])
    [0, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900]
    >>> u3 = mda.Universe(GRO, XTC, XTC, all_coordinates=False, continuous=True)
    >>> u3.trajectory
    <ChainReader containing adk_oplsaa.xtc with 10 frames of 47681 atoms>
    >>> print([int(ts.time) for ts in u3.trajectory])
    [0, 100, 200, 300, 400, 500, 600, 700, 800, 900]

TODO: Figure out how to treat `continuous`.

**The following options modify the created Universe:**

* :code:`guess_bonds`: whether to guess connectivity between atoms. (default: False)
* :code:`vdwradii`: a dictionary of :code:`{element: radius}` of van der Waals' radii for use in guessing bonds.
* :code:`transformations`: a function or list of functions for on-the-fly trajectory transformation.
* :code:`in_memory`: whether to load coordinates into memory (default: False)
* :code:`in_memory_step`: only read every nth frame into an in-memory representation. (default: 1)
* :code:`is_anchor`: whether to consider this Universe when unpickling :code:`AtomGroups` (default: True)
* :code:`anchor_name`: the name of this Universe when unpickling :code:`AtomGroups` (default: None, automatically generated)

----------------------------
Constructing from AtomGroups
----------------------------

A new Universe can be created from one or more :code:`AtomGroup` instances with :func:`MDAnalysis.Merge()`. The :code:`AtomGroup` instances can come from different Universes, meaning that this is one way to concatenate selections from different datasets. 

-------------------------
Constructing from scratch
-------------------------

Universe.empty()



Universe properties
===================

Topology properties
===================


Universe methods
===================

