.. -*- coding: utf-8 -*-

============
Trajectories
============

In MDAnalysis, static data is contained in your universe Topology, while dynamic data is drawn from its trajectory at ``Universe.trajectory``. This is typically loaded from a trajectory file and includes information such as:

    * atom coordinates (``Universe.atoms.positions``)
    * box size (``Universe.dimensions``)
    * velocities and forces (if your file format contains the data) (``Universe.atoms.velocities``)

Although these properties look static, they are actually dynamic, and the data contained within can change.
In order to remain memory-efficient, MDAnalysis does not load every frame of your trajectory into memory at once. Instead, a Universe has a state: the particular timestep that it is currently associated with in the trajectory. When the timestep changes, the data in the properties above shifts accordingly.

The typical way to change a timestep is to index it. ``Universe.trajectory`` can be thought of as a list of :class:`~MDAnalysis.coordinates.base.Timestep`\ s, a data structure that holds information for the current time frame. For example, you can query its length.

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD

    u = mda.Universe(PSF, DCD)
    len(u.trajectory)

When a trajectory is first loaded from a file, it is set to the first frame by default.

.. ipython:: python

    print(u.trajectory.ts, u.trajectory.time)

Indexing the trajectory returns the timestep for that frame, and sets the Universe to point to that frame until the timestep next changes.

.. ipython:: python

    u.trajectory[3]

.. ipython:: python

    print('Time of fourth frame', u.trajectory.time)

Many tasks involve applying a function to each frame of a trajectory. For these, you need to iterate through the frames, *even if you don't directly use the timestep*. This is because the act of iterating moves the Universe onto the next frame, changing the dynamic atom coordinates. 

Trajectories can be sliced as below, if you only want to work on a subset of frames.

.. ipython:: python

    protein = u.select_atoms('protein')
    for ts in u.trajectory[:20:4]:
        # the radius of gyration depends on the changing positions
        rad = protein.radius_of_gyration()
        print('frame={}: radgyr={}'.format(ts.frame, rad))
    
Note that after iterating over the trajectory, the frame is always set back to the first frame, even if your loop stopped before the trajectory end.

.. ipython:: python

    u.trajectory.frame

Because MDAnalysis will pull trajectory data directly from the file it is reading from, changes to atom coordinates and box dimensions will not persist once the frame is changed. The only way to make these changes permanent is to load the trajectory into memory, or to write a new trajectory to file for every frame. For example, to set a cubic box size for every frame and write it out to a file::

    with mda.Writer('with_box.trr', 'w', n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            ts.dimensions = [10, 10, 10, 90, 90, 90]
            w.write(u.atoms)
    
    u_with_box = mda.Universe(PSF, 'with_box.trr')


Sometimes you may wish to only transform part of the trajectory, or to not write a file out. In these cases, MDAnalysis supports "on-the-fly" transformations that are performed on a frame when it is read. 



On-the-fly transformations
==========================

An on-the-fly transformation is a function that modifies the data contained in a trajectory :class:`~MDAnalysis.coordinates.base.Timestep`. It is called for each current time step as it is loaded into memory. A transformation function must also return the current :class:`~MDAnalysis.coordinates.base.Timestep`, as transformations are often chained together.

The :mod:`MDAnalysis.transformations` module contains a collection of transformations. For example, :func:`~MDAnalysis.transformations.fit.fit_rot_trans` can perform a mass-weighted alignment on an :class:`~MDAnalysis.core.groups.AtomGroup` to a reference.

.. ipython:: python

    from MDAnalysis.transformations import fit

    protein = u.select_atoms('protein')
    align_transform = fit.fit_rot_trans(protein, protein, weights='mass')
    u.trajectory.add_transformations(align_transform)

Other implemented transformations include functions to :mod:`~MDAnalysis.transformations.translate`, :mod:`~MDAnalysis.transformations.rotate`, :mod:`~MDAnalysis.transformations.fit` an :class:`~MDAnalysis.core.groups.AtomGroup` to a reference, and :mod:`~MDAnalysis.transformations.wrap` or unwrap an :class:`~MDAnalysis.core.groups.AtomGroup` in the unit cell. 

If you need a different transformation, it is easy to implement your own.

----------------------
Custom transformations
----------------------

At its core, a transformation function must only take a :class:`~MDAnalysis.coordinates.base.Timestep` as its input and return the :class:`~MDAnalysis.coordinates.base.Timestep` as the output.

.. ipython:: python

    def up_by_2(ts):
        """Translates atoms up by 2 angstrom"""
        ts.positions += np.array([0.0, 0.0, 0.2])
        return ts
    
    u = mda.Universe(PSF, DCD, transformations=[up_by_2])


If your transformation needs other arguments, you will need to wrap your core transformation with a wrapper function that can accept the other arguments.

.. ipython:: python

    def up_by_x(x):
        """Translates atoms up by x angstrom"""
        def wrapped(ts):
            """Handles the actual Timestep"""
            ts.positions += np.array([0.0, 0.0, float(x)])
            return ts
        return wrapped
    
    # load Universe with transformations that move it up by 7 angstrom
    u = mda.Universe(PSF, DCD, transformations=[up_by_x(5), up_by_x(2)])

    
Alternatively, you can use :func:`functools.partial` to substitute the other arguments.

.. ipython:: python

    import functools

    def up_by_x(ts, x):
        ts.positions += np.array([0.0, 0.0, float(x)])
        return x
    
    up_by_5 = functools.partial(up_by_x, x=5)
    u = mda.Universe(PSF, DCD, transformations=[up_by_5])

Above we have shown that a :class:`~MDAnalysis.core.universe.Universe` can be created with transformations directly. If your transformation depends on something within the :class:`~MDAnalysis.core.universe.Universe` (e.g. it needs to operate on a particular :class:`~MDAnalysis.core.groups.AtomGroup`), then you can load the :class:`~MDAnalysis.core.universe.Universe` and use the :meth:`~MDAnalysis.core.universe.Universe.add_transformations` method to add transformations.

You can only add transformations *once*, so add your entire workflow at the same time. The below code transforms the trajectory so that the first 100 residues are translated up by 10 angstrom, and the remaining residues are translated down 10 angstrom.

.. ipython:: python

    def ag_up_by_x(ag, x):
        def wrapped(ts):
            ag.positions += np.array([0.0, 0.0, float(x)])
            return ts
        return wrapped

    u = mda.Universe(PSF, DCD)
    res_to_100 = u.residues[:100].atoms
    res_after_100 = u.residues[100:].atoms

    workflow = [ag_up_by_x(res_to_100, 10),
                ag_up_by_x(res_after_100, -10)]
    u.trajectory.add_transformations(*workflow)
    