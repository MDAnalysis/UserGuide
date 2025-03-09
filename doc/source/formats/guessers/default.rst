.. _default-guesser:

==============
DefaultGuesser
==============


.. warning::

    The default guesser has been created to reproduce pre-v2.8.0 MDAnalysis guessing behaviour as much as possible. **Default behaviours will change in MDAnalysis version 3**, as detailed in deprecation warnings.


The :class:`~MDAnalysis.guesser.default_guesser.DefaultGuesser` is the default guessing context for MDAnalysis. For historical reasons the ``DefaultGuesser`` largely works with biological conventions; for example, an atom named CA will be assumed to be a carbon rather than a calcium atom.


Attributes guessed
==================

The topology attributes guessed by the default guesser are listed below, as are their dependencies and broad assumptions.
Please see the `Guesser API documentation`_ for more details.

.. _`Guesser API documentation`: https://docs.mdanalysis.org/stable/documentation_pages/guesser_modules/default_guesser.html

.. _default-guesser-types:

------------------
Elements and types
------------------

The default guesser guesses atom ``element``\ s and ``type``\ s using the same pathway; when atom ``type``\ s are guessed, they represent the atom ``element``. Atom elements are guessed from the atom name. The default guesser follows biological naming conventions, where atoms named "CA" are much more likely to represent an alpha-carbon than a calcium atom. This guesser is still relatively fragile for non-traditionally biological atom names.

The :meth:`~MDAnalysis.guesser.default_guesser.DefaultGuesser.guess_atom_element` method is used to guess atom elements or types following a process by which numbers, symbols, and some letters are stripped from the atom name and checked against a look-up table, as detailed in the `Guesser API documentation`_. With this method, for example, "AO5*" would be guessed as "O", and "3hg2" as "H".


------
Masses
------

Masses are guessed by using a look-up table to get masses from the atom's ``element`` attribute. If ``element``\ s are not available, the atom's ``type`` is used in place of the element. If the ``type`` is not available, that is
:ref:`guessed first <default-guesser-types>`.


.. warning::

    When an atom mass cannot be guessed from the atom ``type`` or ``name``, the atom is currently assigned a mass of 0.0.

    Masses are guessed atom-by-atom, so even if most atoms have been guessed correctly, it is possible that some have been given masses of 0. It is important to check for non-zero masses before using methods that rely on them, such as :meth:`AtomGroup.center_of_mass`.


.. important::

    `np.nan` will be used as a default or "missing" value
    in place of 0.0 for atom masses in version 3.0 of MDAnalysis.


-------------
Aromaticities
-------------

These are guessed using the :ref:`RDKit <RDKit-format>` converter by using the ``GetIsAromatic`` method.

.. note:: 
   RDKit needs to have been installed for aromaticity guessing to be available.
   RDKit is always installed when MDAnalysis was installed with conda-forge packages
   but this may not be the case when using other installation paths.

-----------------------------------
Bonds, Angles, Dihedrals, Impropers
-----------------------------------

MDAnalysis can guess if bonds exist between two atoms, based on the distance between them. A bond is created if the 2 atoms are within

.. math::

    d < f \cdot (R_1 + R_2)

of each other, where :math:`R_1` and :math:`R_2` are the van der Waals radii
of the atoms and :math:`f` is an ad-hoc *fudge factor*. This is
the `same algorithm that VMD uses`_.

.. note::

    Previously, guessing bonds would also guess angles and dihedrals. This is no longer the case. Angles, dihedrals, and impropers are not guessed by the default guesser unless
    explicitly requested by the user.


Angles can be guessed from the bond connectivity. MDAnalysis assumes that if atoms 1 & 2 are bonded, and 2 & 3 are bonded, then (1,2,3) must be an angle.

::

   1
    \
     2 -- 3

Dihedral angles and improper dihedrals can both be guessed from angles. Proper dihedrals are guessed by assuming that if (1,2,3) is an angle, and 3 & 4 are bonded, then (1,2,3,4) must be a dihedral.

::

   1        4
    \      /
     2 -- 3

Likewise, if (1,2,3) is an angle, and 2 & 4 are bonded, then (2, 1, 3, 4) must be an improper dihedral (i.e. the improper dihedral is the angle between the planes formed by (1, 2, 3) and (1, 3, 4))

::

   1
    \
     2 -- 3
    /
   4


.. _`same algorithm that VMD uses`:
    http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html
