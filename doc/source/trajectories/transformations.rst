.. -*- coding: utf-8 -*-
.. _transformations:

On-the-fly transformations
==========================


----------------------
Custom transformations
----------------------

At its core, a transformation function must only take a :class:`~MDAnalysis.coordinates.base.Timestep` as its input and return the :class:`~MDAnalysis.coordinates.base.Timestep` as the output.

.. ipython:: python

    import numpy as np
    
    def up_by_2(ts):
        """Translates atoms up by 2 angstrom"""
        ts.positions += np.array([0.0, 0.0, 0.2])
        return ts
    
    u = mda.Universe(TPR, XTC, transformations=[up_by_2])


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
    u = mda.Universe(TPR, XTC, transformations=[up_by_x(5), up_by_x(2)])

    
Alternatively, you can use :func:`functools.partial` to substitute the other arguments.

.. ipython:: python

    import functools

    def up_by_x(ts, x):
        ts.positions += np.array([0.0, 0.0, float(x)])
        return ts
    
    up_by_5 = functools.partial(up_by_x, x=5)
    u = mda.Universe(TPR, XTC, transformations=[up_by_5])

On-the-fly transformation functions can be applied to any property of a Timestep, not just the atom positions. For example, to give each frame of a trajectory a box:

.. ipython:: python
    
    def set_box(ts):
        # creates box of length 10 on x-axis, 20 on y-axis, 30 on z-axis
        # angles are all 90 degrees
        ts.dimensions = [10, 20, 30, 90, 90, 90]
        return ts
    
    u = mda.Universe(TPR, XTC, transformations=[set_box])

