.. -*- coding: utf-8 -*-
.. _IMD-format:

==================================================================
IMD (Data streamed via Interactive Molecular Dynamics protocol v3)
==================================================================

.. include:: classes/IMD.txt

IMDv2 and IMDv3 enable real-time streaming of simulation data between molecular dynamics engines and receiving clients. The :class:`~MDAnalysis.coordinates.IMD.IMDReader` implements the IMDv3 protocol.

.. note::
   MDAnalysis supports **IMDv3 only**, which provides continuous, gap-free streaming and is implemented in modern versions of GROMACS, LAMMPS, and NAMD. IMDv2, while widely available, was designed primarily for visualization and doesn't enforce a consistent number of integration steps between transmitted frames

What is Streaming?
==================

Streaming involves processing data in real-time as it is generated, rather than storing it for later analysis. In molecular dynamics, this means sending simulation data to a client on-the-fly while the simulation is running, without writing large trajectory files to disk.

In IMDv3, this is achieved through a TCP/IP socket connection between the simulation engine and receiving client, transmitting coordinates, velocities, forces, energies, and timing information.

MDAnalysis's :class:`~MDAnalysis.coordinates.IMD.IMDReader` uses the `imdclient <https://imdclient.readthedocs.io/>`_ package and provides a familiar interface for reading streaming data, similar to other trajectory readers in MDAnalysis.

When to Use Streaming?
======================

Streaming analysis is particularly valuable for:

**Long-running simulations**
    Early detection of problems (crashes, artifacts, equilibration issues) can save computational resources.

**Adaptive sampling workflows**
    Real-time analysis can guide simulation parameters or trigger enhanced sampling methods.

**Interactive research**
    Immediate feedback allows researchers to make informed decisions about continuing, modifying, or terminating simulations.

**Storage-constrained environments**
    Analyze data as it's generated without storing large trajectory files.

Installation and Setup
======================

Required Dependencies
---------------------

The IMDReader requires the ``imdclient`` package:

.. code-block:: bash

   pip install imdclient

.. note::
   MDAnalysis requires ``imdclient >= 0.2.2`` for its current implementation.

MD Engine Configuration
-----------------------

We provide example configurations below for enabling IMDv3 streaming in popular MD engines.

**GROMACS**

Add IMD settings to your ``.mdp`` file:

.. code-block:: text

    ; IMD settings
    IMD-group        = System
    IMD-version      = 3
    IMD-nst          = 1
    IMD-time         = No
    IMD-coords       = Yes
    IMD-vels         = No
    IMD-forces       = No
    IMD-box          = No
    IMD-unwrap       = No
    IMD-energies     = No



Run with IMD enabled:

.. code-block:: bash

    gmx mdrun -v -nt 4 -imdwait -imdport 8889

**LAMMPS**

Use the IMD fix in your input script:

.. code-block:: text

    # IMD setup
    fix ID group-ID imd <port> trate <frequency> version 3 unwrap <on/off> time <on/off> box <on/off> coordinates <on/off> velocities <on/off> forces <on/off>

Run your LAMMPS simulation as usual.

**NAMD**

Add IMD configuration to your NAMD configuration file:

.. code-block:: text

    # IMD Settings
    IMDon yes
    IMDport <port, must be the same port used for the client>
    IMDwait <yes/no>
    IMDfreq <frequency of sending data to the client>

    IMDsendPositions <yes/no>
    IMDsendEnergies <yes/no>
    IMDsendTime <yes/no>
    IMDsendBoxDimensions <yes/no>
    IMDsendVelocities <yes/no>
    IMDsendForces <yes/no>
    IMDwrapPositions <yes/no>

Run your NAMD simulation as usual.

.. seealso::
   For detailed engine-specific setup instructions, see the `imdclient simulation engine documentation <https://imdclient.readthedocs.io/en/latest/usage.html>`_.

Basic Usage
===========

Connecting to a Running Simulation
-----------------------------------

Once your simulation is running with IMD enabled:

.. code-block:: python

   import MDAnalysis as mda
   
   # Connect to the simulation
   u = mda.Universe("topol.tpr", "imd://localhost:8889", buffer_size=10*1024*1024)
   
   # Streaming analysis loop
   for ts in u.trajectory:
       print(f"Time: {ts.time:.2f} ps, Step: {ts.data.get('step', 'N/A')}")
       
       # Your analysis code here
       selected_atoms = u.select_atoms("protein and name CA")
       center_of_mass = selected_atoms.center_of_mass()
       print(f"Protein COM: {center_of_mass}")
       
       # Optional: break on some condition
       if ts.time > 1000:  # Stop after 1000 ps
           break

Real-time Quality Control
-------------------------

Monitor simulation health in real-time:

.. code-block:: python

   import MDAnalysis as mda
   import numpy as np
   
   u = mda.Universe("system.tpr", "imd://localhost:8889")
   
   previous_positions = None
   
   for ts in u.trajectory:
       current_positions = u.atoms.positions.copy()
       
       # Check for simulation artifacts
       if previous_positions is not None:
           displacement = np.linalg.norm(current_positions - previous_positions, axis=1)
           max_displacement = np.max(displacement)
           
           if max_displacement > 10.0:  # Atoms moved > 10 Å in one step
               print(f"WARNING: Large displacement detected at {ts.time} ps: {max_displacement:.2f} Å")
       
       # Monitor energies if available
       if 'potential' in ts.data:
           print(f"Potential energy: {ts.data['potential']:.2f}")
       
       previous_positions = current_positions

Advanced Features
=================

Buffer Management
-----------------

For compute-intensive analysis, increase the buffer size to reduce communication overhead:

.. code-block:: python

   # Larger buffer for better performance
   u = mda.Universe("topol.tpr", "imd://localhost:8889", buffer_size=50*1024*1024)

Connection Management
---------------------

Always ensure proper cleanup, especially in interactive environments like Jupyter notebooks and other interactive environments:

.. code-block:: python

   import MDAnalysis as mda
   
   try:
       u = mda.Universe("topol.tpr", "imd://localhost:8889")
       
       for ts in u.trajectory:
           # Your analysis here
           pass
           
   except Exception as e:
       print(f"Error during streaming: {e}")
   finally:
       # Always close the connection
       u.trajectory.close()

Available Data
--------------

The IMDReader provides access to additional simulation data through ``ts.data``:

* ``dt``: Time step size in picoseconds
* ``step``: Current simulation step number  
* Energy terms: ``potential``, ``total``, etc. (engine-dependent)

.. code-block:: python

   for ts in u.trajectory:
       print(f"Step {ts.data.get('step')}: dt={ts.data.get('dt')} ps")
       
       # Available energy terms vary by MD engine
       for key, value in ts.data.items():
           if key not in ['dt', 'step']:
               print(f"  {key}: {value}")

Multiple Client Connections
===========================

The ability to connect multiple clients to the same IMD port depends on the MD engine implementation:

* **GROMACS**: Typically supports single client connections
* **LAMMPS**: May support multiple clients (version-dependent)
* **NAMD**: Supports multiple clients

.. important::
   Even when multiple connections are supported, each receives an independent data stream. Different clients may receive different data depending on the engine configuration.

Integration with MDAnalysis Tools
=================================

Most MDAnalysis analysis classes work with streaming data, but some limitations apply:

**Compatible Analysis**

.. code-block:: python

   from MDAnalysis.analysis import distances, contacts
   
   u = mda.Universe("system.tpr", "imd://localhost:8889")
   
   for ts in u.trajectory:
       # Distance calculations work normally
       protein = u.select_atoms("protein")
       rg = protein.radius_of_gyration()
       
       # Contact analysis
       selection1 = u.select_atoms("resid 1-10")
       selection2 = u.select_atoms("resid 50-60")
       dist_array = distances.distance_array(selection1.positions, selection2.positions)

**Limitations with Streaming**

Some analysis methods require the complete trajectory and won't work with streaming:

.. code-block:: python

   # These will NOT work with streaming:
   # - trajectory.timeseries()
   # - Most analysis classes that need multiple passes
   # - Random frame access (trajectory[10])
   # - Backward iteration

Important Limitations
=====================

Streaming analysis has fundamental constraints due to its real-time nature:

**Data Access Limitations**

* **No random access**: Cannot jump to arbitrary frames or seek backwards
* **Forward-only**: Can only iterate through frames as they arrive
* **Single-use**: Cannot restart iteration once the stream is consumed
* **No trajectory length**: Total frame count unknown until simulation ends
* **No independent copies**: Cannot create multiple reader instances for the same stream

**Analysis Constraints**

* **No timeseries methods**: Cannot use ``trajectory.timeseries()``
* **No bulk operations**: Cannot extract all data at once
* **Limited multiprocessing**: Cannot split across processes
* **Single client**: Only one reader per IMD stream (engine-dependent)

**Practical Considerations**

.. code-block:: python

   # This WILL work - forward iteration
   for ts in u.trajectory:
       analysis_data.append(calculate_something(ts))
   
   # This will NOT work - random access
   ts = u.trajectory[10]  # TypeError
   
   # This will NOT work - backwards iteration  
   for ts in u.trajectory[::-1]:  # ValueError
       pass
   
   # This will NOT work - restarting iteration
   for ts in u.trajectory:
       break
   for ts in u.trajectory:  # Won't start from beginning
       pass

See Also
========

* :class:`~MDAnalysis.coordinates.IMD.IMDReader` - Technical API documentation
* :class:`~MDAnalysis.coordinates.base.StreamReaderBase` - Base class for streaming readers
* `imdclient documentation <https://imdclient.readthedocs.io/>`_ - Complete imdclient package documentation
* `IMDv3 protocol specification <https://imdclient.readthedocs.io/en/latest/protocol_v3.html>`_ - Technical protocol details