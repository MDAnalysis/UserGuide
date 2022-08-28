.. -*- coding: utf-8 -*-
.. _universe:


Universe
========

    *If you wish to make an apple pie from scratch, you must first invent the universe.*

    -- Carl Sagan, Cosmos

MDAnalysis is structured around two fundamental classes: the :class:`~MDAnalysis.core.universe.Universe` and the :class:`~MDAnalysis.core.groups.AtomGroup`. Almost all code in MDAnalysis begins with :class:`~MDAnalysis.core.universe.Universe`, which contains all the information describing a molecular dynamics system.

It has two key properties:

* :attr:`~MDAnalysis.core.universe.Universe.atoms`: an :class:`~MDAnalysis.core.groups.AtomGroup` of the system's atoms, providing access to important analysis methods (described below)
* :attr:`~MDAnalysis.core.universe.Universe.trajectory`: the currently loaded trajectory reader

A :class:`~MDAnalysis.core.universe.Universe` ties the static information from the "topology" (e.g. atom identities) to dynamically updating information from the "trajectory" (e.g. coordinates). A key feature of MDAnalysis is that an entire trajectory is not loaded into memory (unless the user explicitly does so with :class:`~MDAnalysis.coordinates.memory.MemoryReader`). Instead, the :attr:`~MDAnalysis.core.universe.Universe.trajectory` attribute provides a view on a specific frame of the trajectory. This allows the analysis of arbitrarily long trajectories without a significant impact on memory.

-------------------
Creating a Universe
-------------------

.. _universe-loading:

Loading from files
------------------

A Universe is typically created from a "topology" file, with optional "trajectory" file/s. Trajectory files must have the coordinates in the same order as atoms in the topology. See :ref:`Formats <formats>` for the topology and trajectory formats supported by MDAnalysis, and how to load each specific format. ::

    u = Universe(topology, trajectory)
    u = Universe(pdbfile)                       # read atoms and coordinates from PDB or GRO
    u = Universe(topology, [traj1, traj2, ...]) # read from a list of trajectories
    u = Universe(topology, traj1, traj2, ...)   # read from multiple trajectories

The line between topology and trajectory files is quite blurry. For example, a PDB or GRO file is considered both a topology and a trajectory file. The difference is that a **topology file** provides static information, such as atom identities (name, mass, etc.), charges, and bond connectivity. A **trajectory file** provides dynamic information, such as coordinates, velocities, forces, and box dimensions.

If only a single file is provided, MDAnalysis tries to read both topology and trajectory information from it. When multiple trajectory files are provided, coordinates are loaded in the order given.

The default arguments should create a Universe suited for most analysis applications. However, the :class:`~MDAnalysis.core.universe.Universe` constructor also takes optional arguments.


**The following options specify how to treat the input:**

* :code:`format`: the file format of the trajectory file/s. (default: None, formats are guessed)
* :code:`topology_format`: the file format of the topology file. (default: None, formats are guessed)
* :code:`all_coordinates`: whether to read coordinate information from the first file (default: False. Ignored when only one file is provided)
* :code:`continuous`: whether to give multiple trajectory files continuous time steps. This is currently only supported for XTC/TRR trajectories with a GRO/TPR topology, following the behaviour of `gmx trjcat <http://manual.gromacs.org/documentation/2018/onlinehelp/gmx-trjcat.html>`_ (default: False.)

.. ipython:: python
    :okwarning:

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PDB, GRO, XTC
    u1 = mda.Universe(GRO, XTC, XTC, all_coordinates=True)
    u1.trajectory
    u2 = mda.Universe(GRO, XTC, XTC, all_coordinates=False, continuous=False)
    print([int(ts.time) for ts in u2.trajectory])



**The following options modify the created Universe:**

* :code:`guess_bonds`: whether to guess connectivity between atoms. (default: False)
* :code:`vdwradii`: a dictionary of :code:`{element: radius}` of van der Waals' radii for use in guessing bonds.
* :code:`transformations`: a function or list of functions for on-the-fly trajectory transformation.
* :code:`in_memory`: whether to load coordinates into memory (default: False)
* :code:`in_memory_step`: only read every nth frame into an in-memory representation. (default: 1)
* :code:`is_anchor`: whether to consider this Universe when unpickling :class:`~MDAnalysis.core.groups.AtomGroup`\ s (default: True)
* :code:`anchor_name`: the name of this Universe when unpickling :class:`~MDAnalysis.core.groups.AtomGroup`\ s (default: None, automatically generated)

.. _universe-kwargs:

You can also pass in keywords for parsing the topology or coordinates. For example, many file formats do not specify the timestep for their trajectory. In these cases, MDAnalysis assumes that the default timestep is 1 ps. If this is incorrect, you can pass in a ``dt`` argument to modify the timestep. **This does not modify timesteps for formats that include time information.**

.. ipython:: python
    :okwarning:

    from MDAnalysis.tests.datafiles import PRM, TRJ
    default_timestep = mda.Universe(PRM, TRJ)
    default_timestep.trajectory.dt
    user_timestep = mda.Universe(PRM, TRJ, dt=5)  # ps
    user_timestep.trajectory.dt




Constructing from AtomGroups
----------------------------

A new Universe can be created from one or more :class:`~MDAnalysis.core.groups.AtomGroup` instances with :func:`~MDAnalysis.core.universe.Merge()`. The :class:`~MDAnalysis.core.groups.AtomGroup` instances can come from different Universes, meaning that this is one way to concatenate selections from different datasets.

For example, to combine a protein, ligand, and solvent from separate PDB files:

.. code-block:: python

    u1 = mda.Universe("protein.pdb")
    u2 = mda.Universe("ligand.pdb")
    u3 = mda.Universe("solvent.pdb")
    u = Merge(u1.select_atoms("protein"), u2.atoms, u3.atoms)
    u.atoms.write("system.pdb")


Constructing from scratch
-------------------------

A Universe can be constructed from scratch with :meth:`Universe.empty <MDAnalysis.core.universe.Universe.empty>`. There are three stages to this process:

    #. Create the blank Universe with specified number of atoms. If coordinates, set :code:`trajectory=True`.
    #. Add topology attributes such as atom names.
    #. (Optional) Load coordinates.

For example, to construct a universe with 6 atoms in 2 residues:

.. ipython:: python
    :okwarning:

    u = mda.Universe.empty(6, 2, atom_resindex=[0, 0, 0, 1, 1, 1], trajectory=True)
    u.add_TopologyAttr('masses')
    coordinates = np.empty((1000,  # number of frames
                            u.atoms.n_atoms,
                            3))
    u.load_new(coordinates, order='fac')


`See this notebook tutorial for more information. <examples/constructing_universe.ipynb>`_

.. _guessing-topology-attributes:

----------------------------
Guessing topology attributes
----------------------------

MDAnalysis has a guesser library that hold various guesser classes. Each guesser class is tailored to be context-specific. For example, PDBGuesser is specific for guessing attributes for PDB file format. See :ref:`guessing` for more details about the available context-aware guessers.
The Universe has the ability to guess an attribute within a specific context at the universe creation or by using :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttributes` API.
For example, to guess ``element`` attribute for a PDB file by either of two ways:

.. ipython:: python
    :okwarning:

    u = mda.Universe(PDB, context='PDB', to_guess=['elements'])

or by using the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttributes` API

.. ipython:: python
    :okwarning:

    u = mda.Universe(PDB)
    u.guess_TopologyAttributes(context='PDB', to_guess=['elements'])

**The following options modify how to guess attribute(s):**

* :code:`context`: the context of the guesser to be used in guessing the attribute. You can pass either a string representing the context (see :ref:`guessing` for more detail about available guessers and their context), or as an object of a guesser class. The default value of the context is :code:`default`, which corresponds to a generic :class:`~MDAnalysis.guesser.default_guesser.DefaultGuesser`, that is not specific to any context. You can pass a context once, and whenever you call :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttributes` again without passing context, it will assume that  you still using the same context.
* :code:`to_guess`: a list of the desired attributes to be guessed. This has to be the plural name of the attributes (masses not mass).
* :code:`**kwargs`: to pass any supplemental data to the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttributes` API that can be useful in guessing some attributes (eg. passing vdwradii for bond guessing).

For now, MDAnalysis automatically guess :code:`types` and :code:`masses` for nearly all topology formats at the universe creation level. This is done using the :class:`~MDAnalysis.guesser.default_guesser.DefaultGuesser`, which maybe not accurate for all contexts.
You can easily override this as mentioned above, by using context specific guesser if there is one exists for your data.

.. _universe-properties:

-------------------------------
Universe properties and methods
-------------------------------

A Universe holds master groups of atoms and topology objects:

    * :attr:`atoms`: all Atoms in the system, in an :ref:`AtomGroup <atomgroup>`.
    * :attr:`residues`: all Residues in the system
    * :attr:`segments`: all Segments in the system
    * :attr:`bonds`: all bond TopologyObjects in the system
    * :attr:`angles`: all angle TopologyObjects in the system
    * :attr:`dihedrals`: all dihedral TopologyObjects in the system
    * :attr:`impropers`: all improper TopologyObjects in the system

:ref:`Residues and Segments <residues-and-segments>` are chemically meaningful groups of Atoms.

Modifying a topology is typically done through the :class:`~MDAnalysis.core.universe.Universe`, which contains several methods for adding properties:

    * :meth:`~MDAnalysis.core.universe.Universe.add_TopologyAttr`
    * :meth:`~MDAnalysis.core.universe.Universe.add_Residue`
    * :meth:`~MDAnalysis.core.universe.Universe.add_Segment`

See :ref:`topology-attributes` for more information on which topology attributes can be added, and `<examples/constructing_universe.ipynb>`_ for examples on adding attributes and Segments.
