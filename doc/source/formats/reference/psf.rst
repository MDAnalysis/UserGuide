.. -*- coding: utf-8 -*-
.. _PSF-label:

===================================================
PSF (CHARMM, NAMD, or XPLOR protein structure file)
===================================================

.. include:: classes/PSF.txt

A `protein structure file (PSF) <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html>`_ contains topology information for CHARMM, NAMD, and XPLOR. The MDAnalysis PSFParser only reads information about atoms, bonds, angles, dihedrals, and impropers. While PSF files can include information on hydrogen-bond donor and acceptor groups, MDAnalysis does not read these in.

.. important:: **Atom ids**

    Although PSF files index atoms from 1 in the file, the MDAnalysis PSFParser subtracts 1 to create atom :code:`ids`. This means that if your atom is numbered 3 in your PSF file, it will have an :code:`Atom.id` of 2 in MDAnalysis.

    Atom :code:`indices` are MDAnalysis derived and always index from 0, no matter the file type.