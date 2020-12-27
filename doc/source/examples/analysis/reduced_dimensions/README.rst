.. -*- coding: utf-8 -*-

===================
Dimension reduction
===================

A molecular dynamics trajectory with :math:`N` atoms can be considered
 a path through :math:`3N`-dimensional molecular configuration space.
 It remains difficult to extract important dynamics or compare trajectory
 similarity from such a high-dimensional space. However, collective motions
 and physically relevant states can often be effectively described with
 low-dimensional representations of the conformational space explored over
 the trajectory. MDAnalysis implements two methods for dimensionality
 reduction. 


For computing similarity, see the tutorials in :ref:`trajectory-similarity`.

.. toctree::
   :maxdepth: 1

   /examples/analysis/reduced_dimensions/pca
   /examples/analysis/reduced_dimensions/diffusion_map
