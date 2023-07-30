.. -*- coding: utf-8 -*-
.. _TNG-format:

========================================
TNG (Trajectory Next Generation)
========================================

.. include:: classes/TNG.txt

TNG (The Next Generation) is a highly flexible and high performance trajectory file format for molecular simulations, particularly designed for GROMACS. It includes support for storing both reduced precision coordinates and a lossless, full precision version simultaneously. This enables high quality scientific analysis alongside long term storage in a format that retains full scientific integrity. TNG supports an arbitrary number of data blocks of various types to be added to the file.

Reading in
==========

MDAnalysis supports reading of TNG files. The TNGReader uses the GROMACS tng library for reading the binary TNG format.

The TNG format includes high level C-style API functions for increased ease of use.
