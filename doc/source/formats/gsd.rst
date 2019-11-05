.. -*- coding: utf-8 -*-
.. _GSD-label:

====================
GSD (HOOMD GSD file)
====================

The HOOMD schema GSD file format can contain both topology and trajectory information (output of `HOOMD-blue`_). 

.. important:: 

    The GSD format was developed to support changing numbers of particles, particle types, particle identities and topologies. However, MDAnalysis currently does not support changing topologies. Therefore, the MDAnalysis reader should only be used for trajectories that keep the particles and topologies fixed. 
    
    A user will only get an error if the number of particles changes from the first time step. MDAnalysis does not currently check for changes in the particle identity or topology, and it does not update these over the trajectory.

.. _HOOMD-blue : http://codeblue.umich.edu/hoomd-blue/index.html

.. note:: **Residue resnames**

    Unlike other formats, MDAnalysis treats residue :code:`resnames` from GSD files as integers. These are identical to :code:`resids` and :code:`resnums`. 