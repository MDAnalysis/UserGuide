.. -*- coding: utf-8 -*-
.. _XTC-label:

========================================
XTC (GROMACS compressed trajectory file)
========================================

.. include:: classes/XTC.txt

The GROMACS XTC trajectory compresses data with reduced precision (3 decimal palces by default). MDAnalysis can only read coordinates from these files. See :ref:`TRR-label` for uncompressed files that provide velocity and force information.


Reading in
==========

MDAnalysis uses XDR based readers for GROMACS formats, which store offsets on the disk. The offsets are used to enable access to random frames efficiently. These offsets will be generated automatically the first time the trajectory is opened, and offsets are generally stored in hidden ``*_offsets.npz`` files. Occasionally, MDAnalysis fails to `read XDR offsets, resulting in a fatal error`_.

.. _`read XDR offsets, resulting in a fatal error`: https://github.com/MDAnalysis/mdanalysis/issues/1893

Trajectories split across multiple files can be :ref:`read continuously into MDAnalysis <chainreader-label>` with ``continuous=True``, in the style of `gmx trjcat`_.

.. _`gmx trjcat`: http://manual.gromacs.org/documentation/2018/onlinehelp/gmx-trjcat.html