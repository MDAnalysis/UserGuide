.. -*- coding: utf-8 -*-
.. _auxiliary:

===============
Auxiliary files
===============


Auxiliary readers allow you to read in timeseries data accompanying a trajectory that is not stored in the regular trajectory file.

Supported formats
-----------------

   +-------------------------+--------+-----------+---------------------------------------+
   | Reader                  | Format | Extension | Remarks                               |
   |                         |        | (if file) |                                       |
   +=========================+========+===========+=======================================+
   | :class:`.XVGReader`     | XVG    | xvg       | Produced by Gromacs during simulation |
   |                         |        |           | or analysis.                          |
   |                         |        | (default) | Reads full file on initialisation.    |
   +-------------------------+--------+-----------+---------------------------------------+
   | :class:`.XVGFileReader` | XVG-F  | xvg       | Alternate xvg file reader, reading    |
   |                         |        |           | each step from the file in turn for a |
   |                         |        |           | lower memory footprint.               |
   |                         |        |           | :class:`.XVGReader` is the default    |
   |                         |        |           | reader for .xvg files.                |
   +-------------------------+--------+-----------+---------------------------------------+
   | :class:`.EDRReader`     | EDR    | edr       | Produced by Gromacs during simulation |
   |                         |        |           | or analysis.                          |
   |                         |        |           | Reads full file on initialisation.    |
   +-------------------------+--------+-----------+---------------------------------------+



XVG Files
---------

Reading data directly
---------------------

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import XVG_BZ2  # cobrotoxin protein forces
    aux = mda.auxiliary.core.auxreader(XVG_BZ2)
    aux

In stand-alone use, an auxiliary reader allows you to iterate over each step in a set of auxiliary data.

.. ipython:: python

    for step in aux:
        print(step.data)

Use slicing to skip steps.

.. ipython:: python

    for step in aux[1:2]:
        print(step.time)

The :func:`~MDAnalysis.auxiliary.core.auxreader` function uses the :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` to return an appropriate class. This can guess the format either from a filename, '

.. ipython:: python

    mda.auxiliary.core.get_auxreader_for(XVG_BZ2)

or return the reader for a specified format.

.. ipython:: python

    mda.auxiliary.core.get_auxreader_for(format='XVG-F')


Loading data into a Universe
----------------------------

Auxiliary data may be added to a trajectory Reader through the
:meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary` method. Auxiliary data
may be passed in as a AuxReader instance, or directly as e.g. a filename, in
which case :func:`~MDAnalysis.auxiliary.core.get_auxreader_for` is used to
guess an appropriate reader.

.. ipython:: python
    :okwarning:

    from MDAnalysis.tests.datafiles import PDB_xvf, TRR_xvf
    u = mda.Universe(PDB_xvf, TRR_xvf)
    u.trajectory.add_auxiliary('protein_force', XVG_BZ2)
    for ts in u.trajectory:
        print(ts.aux.protein_force)

Passing arguments to auxiliary data
-----------------------------------

For alignment with trajectory data, auxiliary readers provide methods to
assign each auxiliary step to the nearest trajectory timestep, read all steps
assigned to a trajectory timestep and calculate 'representative' value(s) of
the auxiliary data for that timestep.

To set a timestep or ??

'Assignment' of auxiliary steps to trajectory timesteps is determined from the time
of the auxiliary step, ``dt`` of the trajectory and time at the first frame of the
trajectory. If there are no auxiliary steps assigned to a given timestep (or none within
``cutoff``, if set), the representative value(s) are set to ``np.nan``.

Iterating over auxiliary data
-----------------------------

Auxiliary data may not perfectly line up with the trajectory, or have missing data.

.. ipython:: python
    :okwarning:

    from MDAnalysis.tests.datafiles import PDB, TRR
    u_long = mda.Universe(PDB, TRR)
    u_long.trajectory.add_auxiliary('protein_force', XVG_BZ2, dt=200)
    for ts in u_long.trajectory:
        print(ts.time, ts.aux.protein_force[:4])

The trajectory :class:`~MDAnalysis.coordinates.base.ProtoReader` methods
:meth:`~MDAnalysis.coordinates.base.ProtoReader.next_as_aux` and
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_as_aux` allow for movement
through only trajectory timesteps for which auxiliary data is available.

.. ipython:: python

    for ts in u_long.trajectory.iter_as_aux('protein_force'):
        print(ts.time, ts.aux.protein_force[:4])

This may be used to
avoid representative values set to ``np.nan``, particularly when auxiliary data
is less frequent.

Sometimes the auxiliary data is longer than the trajectory.

.. ipython:: python
    :okwarning:

    u_short = mda.Universe(PDB)
    u_short.trajectory.add_auxiliary('protein_force', XVG_BZ2)
    for ts in u_short.trajectory:
        print(ts.time, ts.aux.protein_force)


In order to acess auxiliary values at every individual step, including those
outside the time range of the trajectory,
:meth:`~MDAnalysis.coordinates.base.ProtoReader.iter_auxiliary` allows iteration
over the auxiliary independent of the trajectory.

.. ipython:: python
    :okwarning:

    for step in u_short.trajectory.iter_auxiliary('protein_force'):
        print(step.data)


To iterate over only a certain section of the auxiliary:

.. ipython:: python
    :okwarning:

    for step in u_short.trajectory.iter_auxiliary('protein_force', start=1, step=2):
        # every 2nd step from 1
        print(step.time)

The trajectory remains unchanged, and the auxiliary will be returned to the current
timestep after iteration is complete.

Accessing auxiliary attributes
------------------------------

To check the values of attributes of an added auxiliary, use
:meth:`~MDAnalysis.coordinates.base.ProtoReader.get_aux_attribute`.

.. ipython:: python

    u.trajectory.get_aux_attribute('protein_force', 'dt')


If attributes are settable, they can be changed using
:meth:`~MDAnalysis.coordinates.base.ProtoReader.set_aux_attribute`.

.. ipython:: python

    u.trajectory.set_aux_attribute('protein_force', 'data_selector', [1])


The auxiliary may be renamed using ``set_aux_attribute`` with 'auxname', or more
directly by using :meth:`~MDAnalysis.coordinates.base.ProtoReader.rename_aux`.

.. ipython:: python

    u.trajectory.ts.aux.protein_force
    u.trajectory.rename_aux('protein_force', 'f')
    u.trajectory.ts.aux.f

Recreating auxiliaries
----------------------

To recreate an auxiliary, the set of attributes necessary to replicate it can
first be obtained with :meth:`~MDAnalysis.auxiliary.base.AuxReader.get_description`.
The returned dictionary can then be passed to
:func:`~MDAnalysis.auxiliary.core.auxreader` to load a new copy of the
original auxiliary reader.

.. ipython:: python

    description = aux.get_description()
    list(description.keys())
    del aux
    mda.auxiliary.core.auxreader(**description)

The 'description' of any or all the auxiliaries added to a trajectory can be
obtained using :meth:`~MDAnalysis.coordinates.base.ProtoReader.get_aux_descriptions`.

.. ipython:: python

    descriptions = u.trajectory.get_aux_descriptions(['f'])

To reload, pass the dictionary into :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_auxiliary`.

.. ipython:: python

    u2 = mda.Universe(PDB, TRR)
    for desc in descriptions:
        u2.trajectory.add_auxiliary(**desc)




EDR Files
---------
EDR files are created by GROMACS during simulations and contain additional time-series non-trajectory data on the
system, such as energies, temperature, or pressure. The EDRReader allows direct reading of these binary files
and associating of the data with trajectory time steps just like the XVGReader does. As all functionality of the base AuxReaders (see XVGReader above) also works with
the EDRReader, this section will highlight functionality unique to the EDRReader.


Standalone Usage
----------------
This a test, please remove me, thank you!

.. ipython:: python

    import MDAnalysis as mda
    import MDAnalysisTests as mdat
    from MDAnalysisTests.datafiles import AUX_EDR
    print(mda.__version__)
    print(mdat.__version__)
    !ls /home/docs/checkouts/readthedocs.org/user_builds/mdanalysisuserguide/conda/230/lib/python3.10/site-packages/MDAnalysisTests/data/
    print(AUX_EDR)

The EDRReader is initialised by passing the path to an EDR file to it.

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysisTests.datafiles import AUX_EDR
    aux = mda.auxiliary.EDR.EDRReader(AUX_EDR)

Dozens of terms can be defined in EDR files. A list of available data is conveniently available unter the `terms` attribute.

.. ipython:: python

    aux.terms

To extract data for plotting, the `get_data` method can be used. It can be used to extract a single, multiple, or all terms as follows:

.. ipython:: python

    temp = aux.get_data("Temperature")
    print(temp.keys())
    energies = aux.get_data(["Potential", "Kinetic En."])
    print(energies.keys())
    all_data = aux.get_data()
    print(all_data.keys())

The times of each data point are always part of the returned dictionary to facilicate plotting.

.. ipython:: python

    import matplotlib.pyplot as plt
    plt.plot(temp["Time"], temp["Temperature"])
    plt.show()


Unit Handling
-------------

On object creation, the EDRReader creates a ``unit_dict`` attribute which contains
information on the units of the data stored within. These units are read from the
EDR file automatically and by default converted to MDAnalysis base units where such unit
is defined. The automatic unit conversion can be disabled by setting the ``convert_units`` kwarg to False.

.. ipython:: python

    aux.unit_dict["Box-X"]
    aux_native = mda.auxiliary.EDR.EDRReader(AUX_EDR, convert_units=False)
    aux_native.unit_dict["Box-X"]


Use with Trajectories
---------------------

An arbitrary number of terms to be associated with a trajectory can be specified as a dictionary.
The dictionary is chosen so that the name to be used in MDAnalysis to access the data is mapped to the name in `aux.terms`.

.. ipython:: python

    from MDAnalysisTests.datafiles import AUX_EDR_TPR, AUX_EDR_XTC
    u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
    term_dict = {"temp": "Temperature", "epot": "Potential"}
    u.trajectory.add_auxiliary(term_dict, aux)
    u.trajectory.ts.aux.epot

Adding all data is possible by omitting the dictionary as follows. It is then necessary to
specify that the EDRReader should be passed as the `auxdata` argument as such:

.. ipython:: python

    u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
    u.trajectory.add_auxiliary(auxdata=aux)
    u.trajectory.aux_list


Selecting Trajectory Frames Based on Auxiliary Data
---------------------------------------------------

One use case for the new auxiliary readers is the selection of frames based on auxiliary data.
To select only those frames with a potential energy below a certain threshold, the
following can be used:

.. ipython:: python

    u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
    term_dict = {"epot": "Potential"}
    u.trajectory.add_auxiliary(term_dict, aux)
    selected_frames = np.array([ts.frame for ts in u.trajectory if ts.aux.epot < -524600])

A slice of the trajectory can then be obtained from the list of frames as such:

.. ipython:: python

    trajectory_slice = u.trajectory[selected_frames]
    print(len(u.trajectory))
    print(len(trajectory_slice))


Memory Usage
------------

It is assumed that the EDR file is small enough to be read in full. However,
since one EDRReader instance is created for each term added to a trajectory,
memory usage monitoring was implemented. A warning will be issued if 1 GB of
memory is used by auxiliary data. This warning limit can optionally be changed
by defining the `memory_limit` (in bytes) when adding data to a trajectory. Below,
the memory limit is set to 200 MB.

.. ipython:: python

    u = mda.Universe(AUX_EDR_TPR, AUX_EDR_XTC)
    term_dict = {"temp": "Temperature", "epot": "Potential"}
    u.trajectory.add_auxiliary(term_dict, aux, memory_limit=2e+08)
