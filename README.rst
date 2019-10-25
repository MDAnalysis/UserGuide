================================
  MDAnalysis User Guide
================================

|build| 

|usergroup| |developergroup|


MDAnalysis_ is a Python toolkit to analyze molecular dynamics
trajectories generated by a wide range of popular simulation packages
including DL_Poly, CHARMM, Amber, NAMD, LAMMPS, and Gromacs. (See the
lists of  supported `trajectory formats`_ and `topology formats`_.)

.. code:: python

   import MDAnalysis as mda

   # Load simulation results with a single line
   u = mda.Universe('topol.tpr','traj.trr')

   # Select atoms
   ag = u.select_atoms('name OH')

   # Atom data made available as Numpy arrays
   ag.positions
   ag.velocities
   ag.forces

   # Iterate through trajectories
   for ts in u.trajectory:
       print(ag.center_of_mass())

There are also a number of tutorials_ on the MDAnalysis homepage that explain
how to conduct RMSD calculations, Alignment and more features of MDAnalysis.



.. Footnotes


.. _trajectory formats: https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1
.. _topology formats: https://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats
.. _Issue 87: https://github.com/MDAnalysis/mdanalysis/issues/87
.. _MDAnalysis: https://www.mdanalysis.org
.. _LICENSE: https://github.com/MDAnalysis/mdanalysis/blob/master/LICENSE
.. _`#286`: https://github.com/MDAnalysis/mdanalysis/issues/286
.. _`MDAnalysis.analysis`: https://docs.mdanalysis.org/documentation_pages/analysis_modules.html
.. _`tutorials`: https://www.mdanalysis.org/pages/learning_MDAnalysis/
.. _`guide`: https://github.com/MDAnalysis/mdanalysis/wiki/Guide-for-Developers

.. |usergroup| image:: https://img.shields.io/badge/Google%20Group-Users-lightgrey.svg
   :alt: User Google Group
   :target: http://users.mdanalysis.org

.. |developergroup| image:: https://img.shields.io/badge/Google%20Group-Developers-lightgrey.svg
   :alt: Developer Google Group
   :target: http://developers.mdanalysis.org


.. |build| image:: https://travis-ci.com/MDAnalysis/UserGuide.svg?branch=master
   :alt: Build Status
   :target: https://travis-ci.com/MDAnalysis/UserGuide