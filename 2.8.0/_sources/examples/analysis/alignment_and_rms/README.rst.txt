.. -*- coding: utf-8 -*-

.. _alignment-and-rms:

==========================
Alignments and RMS fitting
==========================

The :mod:`MDAnalysis.analysis.align` and :mod:`MDAnalysis.analysis.rms`
modules contain the functions used for aligning structures,
aligning trajectories, and calculating root mean squared quantities.

Demonstrations of alignment are in `align_structure`_,
`align_trajectory_first`_, and `align_trajectory`_. Another example of
generating an average structure from an alignment is demonstrated in
`rmsf`_. Typically, trajectories need to be aligned for RMSD and
RMSF values to make sense.

.. note::

    These modules use the fast QCP algorithm to calculate the root mean
    square distance (RMSD) between two coordinate sets :cite:`theobald_rapid_2005` and
    the rotation matrix *R* that minimizes the RMSD :cite:`liu_fast_2009`. Please
    cite these references when using these modules.

.. toctree::
   :maxdepth: 1

   /examples/analysis/alignment_and_rms/aligning_structure_to_another
   /examples/analysis/alignment_and_rms/aligning_trajectory
   /examples/analysis/alignment_and_rms/aligning_trajectory_to_frame
   /examples/analysis/alignment_and_rms/rmsd
   /examples/analysis/alignment_and_rms/pairwise_rmsd
   /examples/analysis/alignment_and_rms/rmsf


.. _align_structure: /examples/analysis/alignment_and_rms/aligning_structure_to_another.ipynb
.. _align_trajectory_first: /examples/analysis/alignment_and_rms/aligning_trajectory_to_first_frame.ipynb
.. _align_trajectory: /examples/analysis/alignment_and_rms/aligning_trajectory.ipynb
.. _rmsd: /examples/analysis/alignment_and_rms/rmsd.ipynb
.. _pairwise: /examples/analysis/alignment_and_rms/pairwise_rmsd.ipynb
.. _rmsf: /examples/analysis/alignment_and_rms/rmsf.ipynb
