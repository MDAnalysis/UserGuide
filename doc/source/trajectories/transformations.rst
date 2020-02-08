.. -*- coding: utf-8 -*-
.. _transformations:

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
    