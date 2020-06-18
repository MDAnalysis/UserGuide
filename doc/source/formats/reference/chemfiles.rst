.. -*- coding: utf-8 -*-
.. _CHEMFILES-format:

========================================
CHEMFILES (chemfiles Trajectory or file)
========================================

.. include:: classes/CHEMFILES.txt

The chemfiles_ library provides a C++ implementation of readers and writers for multiple formats. Pass in either a :class:`chemfiles.Trajectory` to be converted into an MDAnalysis Universe, or pass in files to this format with the keyword ``format='CHEMFILES'`` to read the information with the chemfiles implementation. You can also write the MDAnalysis Universe back into a chemfiles object for further work.

.. _chemfiles: https://chemfiles.org/