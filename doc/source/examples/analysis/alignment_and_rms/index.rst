.. -*- coding: utf-8 -*-

.. _alignment_and_rms:

==========================
Alignments and RMS fitting
==========================

The :mod:`MDAnalysis.analysis.align` and :mod:`MDAnalysis.analysis.rms` modules contain the functions used for aligning structures, aligning trajectories, and calculating root mean squared quantities. 

.. note::

    These modules use the fast QCP algorithm to calculate the root mean
    square distance (RMSD) between two coordinate sets [Theobald2005]_ and
    the rotation matrix *R* that minimizes the RMSD [Liu2010]_. Please
    cite these references when using these modules.


.. toctree::
   :maxdepth: 1

   /examples/analysis/alignment_and_rms/aligning_structure_to_another
   /examples/analysis/alignment_and_rms/aligning_structure_to_another.nbconvert
   /examples/analysis/alignment_and_rms/aligning_trajectory_to_first_frame
   /examples/analysis/alignment_and_rms/aligning_trajectory
   /examples/analysis/alignment_and_rms/rmsd
   /examples/analysis/alignment_and_rms/pairwise_rmsd
   /examples/analysis/alignment_and_rms/rmsf