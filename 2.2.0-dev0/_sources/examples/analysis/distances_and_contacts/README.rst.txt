.. -*- coding: utf-8 -*-

======================
Distances and contacts
======================

The :mod:`MDAnalysis.analysis.distances` module provides functions to rapidly compute distances. These largely take in coordinate arrays. 


.. toctree::
   :maxdepth: 1

   /examples/analysis/distances_and_contacts/distances_between_atomgroups
   /examples/analysis/distances_and_contacts/distances_between_selections
   /examples/analysis/distances_and_contacts/distances_within_selection

Residues can be determined to be in contact if atoms from the two residues are within a certain distance. Analysing the fraction of contacts retained by a protein over at trajectory, as compared to the number of contacts in a reference frame or structure, can give insight into folding processes and domain movements. 

:mod:`MDAnalysis.analysis.contacts` contains functions and a class to analyse the fraction of native contacts over a trajectory. 

.. toctree::
   :maxdepth: 1

   /examples/analysis/distances_and_contacts/contacts_native_fraction
   /examples/analysis/distances_and_contacts/contacts_q1q2
   /examples/analysis/distances_and_contacts/contacts_within_cutoff
   /examples/analysis/distances_and_contacts/contacts_custom