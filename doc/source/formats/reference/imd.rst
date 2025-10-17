.. -*- coding: utf-8 -*-
.. _IMD-format:

==================================================================
IMD (Data streamed via Interactive Molecular Dynamics protocol v3)
==================================================================

.. include:: classes/IMD.txt

**Interactive Molecular Dynamics (IMD)** is a data streaming utility that enables transfer of live simulation data (like coordinates, energies, forces) from molecular dynamics engines to receiving clients while the simulation is running. The :class:`~MDAnalysis.coordinates.IMD.IMDReader` is a reader class for accessing IMD live-streamed data by connecting to a live TCP/IP connection.

IMD is particularly useful for long-running simulations, adaptive sampling workflows, real-time quality control, and storage-constrained environments where immediate analysis is needed without writing large trajectory files. More details are provided below.

.. note::
   MDAnalysis supports **IMDv3 only**, which uses the IMDv3 protocol to provide continuous, gap-free streaming and is implemented in modern versions of GROMACS, LAMMPS, and NAMD. IMDv2, while widely available, was designed primarily for visualization and doesn't enforce a consistent number of integration steps between transmitted frames

What is Streaming?
==================

Streaming involves sending and processing data in real-time as it is generated, rather than storing it for later analysis. In molecular dynamics, this means transmitting simulation data to a client on-the-fly while the simulation is running, without writing large trajectory files to disk. The simulation engine acts as a data producer, while the analysis client serves as a receiver that receives this data in real-time.

Interactive Molecular Dynamics (IMD)
------------------------------------

Interactive Molecular Dynamics (IMD) is a specific implementation of streaming for molecular dynamics simulations [#imd2001]_. IMD establishes a TCP/IP socket connection between the simulation engine and receiving client, enabling real-time transmission of simulation data while the simulation is running. It uses a custom protocol that governs the format and type of simulation data that is exchanged.

The previous implementation, IMDv2, was primarily designed for visualization purposes and was built to connect to the `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ molecular visualization program [#imdv2]_. It enabled the transmission of coordinates and energies to VMD for real-time visualization of running simulations. The IMDv3 implementation, however, is designed for more general streaming applications beyond visualization, enabling the transfer of positions, velocities, forces, energies, and timing information [#imdv3]_.

IMDv3 in MD Engines
^^^^^^^^^^^^^^^^^^^^

The IMDv3 protocol is currently implemented in popular MD engines including GROMACS, LAMMPS, and NAMD. Each engine provides configuration options to enable IMDv3 streaming and specify which data types to transmit during simulation execution.

IMDReader in MDAnalysis
^^^^^^^^^^^^^^^^^^^^^^^^

MDAnalysis's :class:`~MDAnalysis.coordinates.IMD.IMDReader` supports the IMDv3 protocol and provides a convenient interface for reading streaming data. It uses the `imdclient <https://imdclient.readthedocs.io/>`_ package which acts as a client/receiver, gathering data from TCP/IP sockets and making it available through MDAnalysis's familiar Universe interface.

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

The IMDReader requires the ``imdclient`` package, which is an optional dependency that must be installed separately from MDAnalysis using either ``pip`` or ``mamba/conda``. For detailed installation instructions, see the `imdclient installation guide <https://imdclient.readthedocs.io/en/stable/getting_started.html#installation>`_.

.. note::
   MDAnalysis requires ``imdclient >= 0.2.2`` for its current implementation.

MD Engine Configuration
-----------------------

To enable IMDv3 streaming, you need to configure your simulation engine with the appropriate IMD settings. The key parameters are:

* **IMD version**: Must be set to 3 for compatibility with MDAnalysis
* **Port**: Network port for the IMD connection (typically 8889)
* **Data transfer rate**: How frequently data is sent (e.g., every 100 steps)
* **Data types**: Which simulation data to stream (coordinates, energies, velocities, etc.)

Below are minimal example configurations for each supported engine. For comprehensive setup instructions and advanced configuration options, see the `imdclient simulation engine documentation <https://imdclient.readthedocs.io/en/latest/usage.html>`_.

GROMACS
^^^^^^^

Add these IMD settings to your ``.mdp`` file:

.. code-block:: text

    ; Required IMD settings
    IMD-version      = 3        ; Use IMDv3 protocol (2/3)
    IMD-group        = System   ; Atom group to stream (System = all atoms)
    IMD-nst          = 100      ; Send data every 100 steps
    
    ; Data streams - specify what to send
    IMD-time         = Yes       ; Stream timing information (Yes/No)
    IMD-energies     = Yes       ; Stream energy data (Yes/No)
    IMD-box          = Yes       ; Stream box dimensions (for PBC analysis) (Yes/No)
    IMD-coords       = Yes       ; Stream coordinates (Yes/No)
    IMD-vels         = Yes       ; Stream velocities (Yes/No)
    IMD-forces       = Yes       ; Stream forces (Yes/No)
    
    ; Coordinate processing
    IMD-unwrap       = Yes       ; Unwrap coordinates across PBC (Yes/No)

Run with IMD enabled:

.. code-block:: bash

    gmx mdrun -v -nt 4 -imdwait -imdport 8889

LAMMPS
^^^^^^

Use the IMD fix in your input script:

.. code-block:: text

    # Complete IMD configuration
    fix imd all imd 8889 trate 100 version 3 unwrap yes time yes box yes coordinates yes velocities yes forces yes
    
    # Parameters explained:
    # 8889: Network port for connection
    # trate 100: Send data every 100 timesteps (transmission rate)
    # version 3: Use IMDv3 protocol (2/3) 
    # unwrap yes: Unwrap coordinates across PBC (yes/no)
    # time yes: Send timing information (yes/no)
    # box yes: Send box dimensions (for PBC analysis) (yes/no)
    # coordinates yes: Send atomic coordinates (yes/no)
    # velocities yes: Send velocity data (yes/no)
    # forces yes: Send force data (yes/no)
    # Note: Energy streaming not supported in LAMMPS IMD fix

Run your LAMMPS simulation as usual.

NAMD
^^^^

Add these IMD settings to your configuration file:

.. code-block:: text

    # Required IMD settings
    IMDon yes               ; Enable IMD functionality
    IMDversion 3            ; Use IMDv3 protocol (2/3)
    IMDport 8889            ; Network port for connection
    IMDwait on              ; Wait for client connection before starting (on/off)
    IMDfreq 100             ; Send data every 100 steps (transmission rate)
    
    # Data streams - specify what to send
    IMDsendTime yes         ; Send timing information (yes/no)
    IMDsendEnergies yes     ; Send energy information (yes/no)
    IMDsendBoxDimensions yes ; Send simulation box data (for PBC analysis) (yes/no)
    IMDsendPositions yes    ; Send coordinates (yes/no)
    IMDsendVelocities yes   ; Send velocity data (yes/no)
    IMDsendForces yes       ; Send force data (yes/no)
    
    # Coordinate processing
    IMDwrapPositions no    ; Don't wrap (i.e unwrap) coordinates into simulation box (yes/no)

Run your NAMD simulation as usual.

Basic Usage
===========

Connecting to a Running Simulation
-----------------------------------

Once your simulation is running with IMD enabled:

.. code-block:: python

   import MDAnalysis as mda
   
   # Connect to the simulation
   u = mda.Universe("topol.tpr", "imd://localhost:8889", buffer_size=10*1024*1024)

   # Select atoms for analysis
   selected_atoms = u.select_atoms("protein and name CA")
   
   # Streaming analysis loop
   for ts in u.trajectory:
       print(f"Time: {ts.time:.2f} ps, Step: {ts.data.get('step', 'N/A')}")
       
       # Your analysis code here
       center_of_mass = selected_atoms.center_of_mass()
       print(f"Protein COM: {center_of_mass}")
       
       # Optional: break on some condition
       if ts.time > 1000:  # Stop after 1000 ps
           break

The ``buffer_size`` parameter (10 MB = 10*1024*1024 bytes in this example) controls the buffer used by `imdclient <https://imdclient.readthedocs.io/>`_ to temporarily store data received from the socket. This buffer accounts for speed differences between the producer (simulation engine) and receiver (analysis code), preventing data loss when analysis is slower than data transmission. A larger buffer is more suitable for systems with many atoms or high transmission frequencies. For more details on buffer management and optimization, see the `imdclient documentation <https://imdclient.readthedocs.io/>`_.

Real-time Quality Control
-------------------------

Monitor simulation health in real-time:

.. code-block:: python

   import MDAnalysis as mda
   from MDAnalysis.lib.distances import calc_bonds
   import numpy as np
   
   # Connect to simulation streaming positions, box dimensions, and energies
   u = mda.Universe("system.tpr", "imd://localhost:8889")
   
   previous_positions = None
   
   for ts in u.trajectory:
       current_positions = u.atoms.positions.copy()
       
       # Check for simulation artifacts using PBC-aware distance calculation
       if previous_positions is not None:
           # Create atom pairs for displacement calculation
           atom_pairs = np.column_stack([np.arange(len(u.atoms)), np.arange(len(u.atoms))])
           
           # Use PBC-aware distance calculation
           displacements = calc_bonds(previous_positions, current_positions, 
                                    atom_pairs, box=ts.dimensions)
           max_displacement = np.max(displacements)
           
           if max_displacement > 10.0:  # Atoms moved > 10 Å in one step
               print(f"WARNING: Large displacement detected at {ts.time} ps: {max_displacement:.2f} Å")
       
       # Monitor energies
       print(f"Potential energy: {ts.data['potential']:.2f}")
       
       previous_positions = current_positions

Advanced Features
=================

Buffer Management
-----------------

The ``buffer_size`` parameter (specified in bytes) controls how much data imdclient can temporarily store while managing speed differences between simulation and analysis. For compute-intensive analysis, increase the buffer size to reduce communication overhead:

.. code-block:: python

   # Larger buffer for demanding scenarios (50 MB = 50*1024*1024 bytes)
   u = mda.Universe("topol.tpr", "imd://localhost:8889", buffer_size=50*1024*1024)

For detailed information about buffer behavior and usage, see the `imdclient buffer management documentation <https://imdclient.readthedocs.io/>`_.

Connection Management
---------------------

Always ensure proper cleanup, especially in interactive environments like Jupyter notebooks and other interactive environments:

.. code-block:: python

   import MDAnalysis as mda
   
   u = None
   error = None
   
   try:
       u = mda.Universe("topol.tpr", "imd://localhost:8889")
       
       for ts in u.trajectory:
           # Your analysis here
           pass
           
   except Exception as e:
       # Log error but don't re-raise yet
       error = e
       print(f"Error during streaming: {e}")

   finally:
       # Always close the connection first
       if u is not None:
           u.trajectory.close()
   
   # Re-raise after cleanup is done
   if error:
       raise error

Available Data
--------------

The IMDReader provides access to additional simulation data through ``ts.data``:

* ``dt``: Time step size in picoseconds
* ``step``: Current simulation step number
* Energy terms: ``potential_energy``, ``total_energy``, etc. (IMD-streamed in NAMD and GROMACS only)

.. code-block:: python

   for ts in u.trajectory:
       print(f"Step {ts.data.get('step')}: dt={ts.data.get('dt')} ps")
       
       # Energy terms available only when IMD-streaming from NAMD and GROMACS only
       for key, value in ts.data.items():
           if key not in ['dt', 'step']:
               print(f"  {key}: {value}")

Integration with MDAnalysis Tools
=================================

Most MDAnalysis analysis classes work with streaming data, but some limitations apply:

Compatible Analysis
^^^^^^^^^^^^^^^^^^^

**What works with streaming:**

* **Single-frame calculations**: Analyses that work on individual timesteps, for example:

  - Within-frame: :meth:`~MDAnalysis.core.groups.AtomGroup.center_of_mass`, :meth:`~MDAnalysis.core.groups.AtomGroup.radius_of_gyration`
  - Within-frame: :func:`~MDAnalysis.analysis.distances.distance_array` for pairwise distances between atom groups
  - Between-frames: :func:`~MDAnalysis.lib.distances.calc_bonds` for displacement calculations comparing consecutive frames
  - Between-frames: Frame-to-frame RMSD calculations using :func:`~MDAnalysis.analysis.rms.rmsd` with stored reference coordinates
  - Within-frame: Real-time monitoring using :meth:`~MDAnalysis.core.groups.AtomGroup.center_of_mass` for quality control checks

* **Accumulative analyses**: Building results incrementally across frames using :class:`~MDAnalysis.analysis.base.AnalysisBase` patterns, for example:

  - :class:`~MDAnalysis.analysis.rdf.InterRDF` - Frame-by-frame radial distribution function calculations
  - :class:`~MDAnalysis.analysis.dihedrals.Dihedral` - Dihedral angle accumulation for conformational analysis  
  - :class:`~MDAnalysis.analysis.lineardensity.LinearDensity` - Density profile building over streaming frames

**What doesn't work:**

* **Multi-pass analyses**: Methods requiring multiple trajectory passes, for example:

  - :class:`~MDAnalysis.analysis.rms.RMSD` - Needs reference structure alignment across all frames / entire trajectory  
  - :class:`~MDAnalysis.analysis.pca.PCA` - Principal component analysis requires full trajectory

* **Global trajectory methods**: Analyses needing simultaneous access to all frames, for example:

  - :meth:`~MDAnalysis.coordinates.base.ProtoReader.timeseries` - Bulk coordinate extraction
  - :class:`~MDAnalysis.analysis.encore.encore` - Ensemble similarity calculations

**Example streaming-compatible analyses:**

.. code-block:: python

   from MDAnalysis.analysis import distances, contacts
   
   u = mda.Universe("system.tpr", "imd://localhost:8889")
   
   # Select atoms once outside the loop (best practice for performance)
   protein = u.select_atoms("protein")
   selection1 = u.select_atoms("resid 1-10")
   selection2 = u.select_atoms("resid 50-60")
   
   for ts in u.trajectory:
       # Distance calculations work normally
       rg = protein.radius_of_gyration()
       
       # Contact analysis
       dist_array = distances.distance_array(selection1.positions, selection2.positions)

Important Limitations
=====================

Streaming analysis has fundamental constraints due to its real-time nature:

Data Access Limitations
^^^^^^^^^^^^^^^^^^^^^^^^

* **No random access**: Cannot jump to arbitrary frames or seek backwards
* **Forward-only**: Can only iterate through frames as they arrive
* **Single-use**: Cannot restart iteration once the stream is consumed
* **No trajectory length**: Total frame count unknown until simulation ends
* **No independent copies**: Cannot create multiple reader instances for the same stream

.. note::
   Multiple client connections to the same IMD port may be possible with some MD engines. For details on engine-specific behavior, see the `imdclient documentation <https://imdclient.readthedocs.io/>`_.

Analysis Constraints
^^^^^^^^^^^^^^^^^^^^

* **No timeseries methods**: Cannot use ``trajectory.timeseries()``
* **No bulk operations**: Cannot extract all data at once
* **Limited multiprocessing**: Cannot split across processes
* **Single client**: Only one reader per IMD stream (engine-dependent)

Practical Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^

**Forward iteration works correctly:**

.. code-block:: python

   # This WILL work - forward iteration
   for ts in u.trajectory:
       analysis_data.append(calculate_something(ts))

**Random frame access will fail:**

.. code-block:: python

   # This will NOT work - random access
   ts = u.trajectory[10]  # ValueError

**Backward iteration will fail:**

.. code-block:: python

   # This will NOT work - backwards iteration  
   for ts in u.trajectory[::-1]:  # ValueError
       pass

**Setting end-frames is not supported:**

.. code-block:: python

   # This will NOT work - cannot set stop index
   for ts in u.trajectory[:10]:  # ValueError
       pass

**Restarting iteration will not work as expected:**

.. code-block:: python

   # This will NOT work - restarting iteration
   for ts in u.trajectory:
       break
   for ts in u.trajectory:  
   # Won't start from beginning but rather continue from where it left off
       pass

See Also
========

* :class:`~MDAnalysis.coordinates.IMD.IMDReader` - Technical API documentation
* :class:`~MDAnalysis.coordinates.base.StreamReaderBase` - Base class for streaming readers
* `imdclient documentation <https://imdclient.readthedocs.io/>`_ - Complete imdclient package documentation
* `IMDv3 protocol specification <https://imdclient.readthedocs.io/en/latest/protocol_v3.html>`_ - Technical protocol details

References
==========

.. [#imd2001] John E. Stone, Justin Gullingsrud, and Klaus Schulten. A system for interactive molecular dynamics simulation. In Proceedings of the 2001 Symposium on Interactive 3D Graphics, I3D '01, 191–194. New York, NY, USA, 2001. Association for Computing Machinery. `<https://dl.acm.org/doi/10.1145/364338.364398>`_

.. [#imdv2] IMDv2 protocol implementation. `VMD Interactive Molecular Dynamics <https://www.ks.uiuc.edu/Research/vmd/imd/>`_. Accessed: 2024.

.. [#imdv3] IMDv3 protocol specification. `imdclient documentation <https://imdclient.readthedocs.io/en/latest/protocol_v3.html>`_. Accessed: 2024.