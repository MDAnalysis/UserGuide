.. -*- coding: utf-8 -*-
.. _guessing:

============================
Guessing Topology Attributes
============================

Since version 2.8.0 MDAnalysis has introduced a new context-dependent guessing API to guess topology attributes that are not read from the file. This allows topology attributes such as masses, charges, and atom types to be guessed from existing information in a context-dependent manner (e.g. biological naming conventions) rather than file formats, as was previously done.

.. list-table:: Supported guesser contexts
   :widths: 25 25 50
   :header-rows: 1

   * - Guesser
     - Context name
     - Topology attributes guessed
   * - :ref:`DefaultGuesser <default-guesser>`
     - "default"
     - elements, types, masses, bonds, angles, dihedrals, impropers, aromaticities


Guessing at Universe creation
=============================

Topology attributes can be guessed at Universe creation by passing in topology attributes to guess to the ``to_guess`` keyword. By default, as of version 2.8.0, the default guesser is used to guess ``types`` and ``masses``.


.. ipython:: python

  import MDAnalysis as mda
  from MDAnalysis.tests.datafiles import PRM12
  u = mda.Universe(PRM12, context="default", to_guess=["types", "masses", "bonds"])
  u.atoms.bonds


In general, guessing at Universe creation works very similarly to guessing using the guess_TopologyAttrs interface documented below. The main difference is that passing guesser-specific keyword arguments such as ``fudge_factor`` and ``vdwradii`` into Universe creation is **now deprecated and will be removed in version 3.0**. Instead, we recommend specifying these arguments through an explicit call to the :ref:`guess_TopologyAttrs <guess-topologyAttrs>` method.

.. _guess-topologyAttrs:

Guessing using the guess_TopologyAttrs interface
================================================

Topology attributes can also be guessed after Universe creation using the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` method. The ``to_guess``, ``force_guess``, and ``context`` keywords are used to specify which attributes to guess, which attributes to forcibly re-guess, and which guesser to use, respectively. These three keywords perform the same way here as they do in Universe creation.

As with Universe creation, the :ref:`DefaultGuesser <default-guesser>` is used as the default ``context``. The following example demonstrates how to guess atom types and masses after Universe creation.

.. ipython:: python

  u = mda.Universe(PRM12, to_guess=[]) # in v2.8.0 masses and types are guessed by default
  u.guess_TopologyAttrs(to_guess=["types", "masses"])
  u.atoms.types


The context can be specified either using a string (e.g. ``"default"``) or an already created Guesser object. It may be convenient to pass in an already-created Guesser object if there are particular keywords you want to use in guessing methods, such as the ``fudge_factor``, ``vdwradii`` or ``lower_bound`` keywords for controlling bond guessing. However, if additional keyword arguments are passed into ``guess_TopologyAttrs``, they will **replace** any existing arguments inside the guesser.

.. ipython:: python

  from MDAnalysis.guesser import DefaultGuesser
  from MDAnalysis.tests.datafiles import CONECT

  u = mda.Universe(CONECT)
  guesser = DefaultGuesser(u, fudge_factor=1.2)
  u.guess_TopologyAttrs(to_guess=["bonds"], context=guesser, fudge_factor=0.5)
  guesser._kwargs["fudge_factor"]


--------------------
Forcibly re-guessing
--------------------

MDAnalysis will preferentially read topology attributes from file instead of re-guessing them, even if the attribute is passed into ``to_guess``. For example, below, the ``types`` attributes reflects the actual atom types in the file.

.. ipython::

  u = mda.Universe(PRM12, to_guess=["types", "masses"])
  u.atoms.types

.. note::

  In cases where the attribute is only present for *some* atoms in the file (e.g. a patchy element column in a PDB), MDAnalysis will only guess the attribute for atoms where it is not present in the file.

To force MDAnalysis to re-guess a TopologyAttr, pass in the attribute to the ``force_guess`` keyword. This will force MDAnalysis to guess the attribute even if it is present in the file.

.. ipython:: python

  u.guess_TopologyAttrs(to_guess=["types"], force_guess=["types"])
  u.atoms.types


------------------------------------
Guessing bonds, angles, and torsions
------------------------------------

Whereas most attributes are guessed at the atom, residue, or segment level, guessing topology objects such as bonds, angles, dihedrals and impropers behaves somewhat differently, and interacts with the ``force_guess`` keyword specially.

Specifically, if these connectivity attributes are guessed, they are by default guessed *additively*. Therefore, if bonds and other objects are guessed twice, the bonds of the second guess are added on. Below, we see the number of bonds increase when guessed again with a looser criteria.

.. ipython:: python

  from MDAnalysis.tests.datafiles import CONECT

  u = mda.Universe(CONECT, to_guess=["bonds"])
  print(len(u.bonds))
  u.guess_TopologyAttrs(to_guess=["bonds"], fudge_factor=1.2) # looser
  print(len(u.bonds))


However, the number of bonds doesn't change when the bonds are guessed again with a stricter criteria -- no new bonds are found.

.. ipython:: python

  u.guess_TopologyAttrs(to_guess=["bonds"], fudge_factor=0.5) # stricter
  print(len(u.bonds))


Moreover, bonds are unique, so if the bonds are guessed again with a same criteria, the guessed bonds doesn't change.

.. ipython:: python

  u.guess_TopologyAttrs(to_guess=["bonds"], fudge_factor=0.5) # same
  print(len(u.bonds))


However, if you want to forcibly overwrite all existing bonds, angles, dihedrals or impropers, you can pass the object to the ``force_guess`` keyword. This will remove all existing objects of that type before guessing. Below, we see the number of bonds has shrunk when guessed with stricter criteria.

.. ipython:: python

  u.guess_TopologyAttrs(to_guess=["bonds"], force_guess=["bonds"], fudge_factor=0.5)
  print(len(u.bonds))


-----------------
Order of guessing
-----------------

The order of the attributes guessed can matter in some cases. For example, bond guessing with the DefaultGuesser relies on looking up the vdW radii of the atoms involved by their atom ``type``. (**Note: this behaviour will likely change to looking up by element in version 3.0**). That means that for file formats where the atom ``type`` is not a valid element, the atom ``type`` must be forcefully re-guessed for bond-guessing to work.

Therefore the following will not work due to the types encoded in the PSF file:

.. ipython:: python
  :okexcept:

  from MDAnalysis.tests.datafiles import PSF, DCD
  u = mda.Universe(PSF, DCD)
  print(u.atoms.types)
  u.guess_TopologyAttrs(to_guess=["bonds"])

However, the snippet below will re-guess the types, and now bond-guessing can work as the elements have vdW radii defined:

.. ipython:: python

  u.guess_TopologyAttrs(to_guess=["types", "bonds"], force_guess=["types"])
  print(u.atoms.types)
