.. -*- coding: utf-8 -*-
.. _guessing:

====================
Guessing
====================

When a Universe is created from a Universe, MDAnalysis can guesses properties that have not been read from the file. Sometimes these properties are available in the file, but are simply not read by MDAnalysis. For example, masses are always guessed.
The :mod:`~MDAnalysis.guesser` module contains different context-specific guessers. This can be forcefield-specific like :class:`~MDAnalysis.guesser.default_guesser.MartiniGuesser`, or format-specific guesser like :class:`~MDAnalysis.guesser.default_guesser.PDBGuesser`.
You can utilize guessers either by initiating an object of it or through the :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttributes` API of the universe to guess various properties to the universe. See :ref:`Guessing topology attributes <guessing-topology-attributes>` for details.

.. _available-guessers:

Available guessers
===================


Here is a list of the currently available context-specific guesser and what attributes they can guess


+--------------------------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------------+
| **guesser**                                | **context**           | **to_guess**                                                                                                    |
+--------------------------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------------+
| :ref:`DefaultGuesser <default-guesser>`    | :code:`default`       | masses, atom types, elements, bonds, angles, dihedrals, improper dihedrals, aromaticities, gasteiger charges    |
+--------------------------------------------+-----------------------+-----------------------------------------------------------------------------------------------------------------+
