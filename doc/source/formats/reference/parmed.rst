.. -*- coding: utf-8 -*-
.. _PARMED-format:

=========================
PARMED (ParmEd Structure)
=========================

.. include:: classes/PARMED.txt

The ParmEd_ library is a general tool for molecular modelling, often used to manipulate system topologies or convert between file formats. You can pass in a :class:`parmed.Structure` to be converted into an MDAnalysis Universe, and convert it back using a :class:`~MDAnalysis.coordinates.ParmEd.ParmEdConverter`. While you can convert Universes created by other means (e.g. by reading files) into a ParmEd structure, MDAnalysis does not read or generate details such the "type" of a bond (i.e. bondtype), or information such as ``ureybradley`` terms. 

.. _ParmEd: https://parmed.github.io/ParmEd/html/index.html