.. currentmodule:: MDAnalysis

=============
API reference
=============

Universe
========

Creating a Universe
-------------------

.. autosummary::
    :toctree: generated/

    Universe
    Merge
    Universe.copy
    Universe.empty
    Universe.load_new


Topology and Group attributes
-----------------------------

.. autosummary::
    :toctree: generated/

    Universe.atoms
    Universe.angles
    Universe.bonds
    Universe.dihedrals
    Universe.impropers
    Universe.residues
    Universe.segments

Other attributes
----------------

.. autosummary::
    :toctree: generated/

    Universe.coord
    Universe.dimensions
    Universe.is_anchor
    Universe.kwargs
    Universe.trajectory


Topology manipulation
---------------------

.. autosummary::
    :toctree: generated/

    Universe.add_Residue
    Universe.add_Segment
    Universe.add_TopologyAttr

Other (rename)
--------------

.. autosummary::
    :toctree: generated/

    Universe.select_atoms
    Universe.transfer_to_memory
    Universe.remove_anchor



AtomGroup
=========

Set methods and group operators
-------------------------------

.. autosummary::
    :toctree: generated/

    AtomGroup.concatenate
    AtomGroup.difference
    AtomGroup.intersection
    AtomGroup.is_strict_subset
    AtomGroup.isdisjoint
    AtomGroup.issubset
    AtomGroup.issuperset
    AtomGroup.subtract
    AtomGroup.symmetric_difference
    AtomGroup.union

Properties
----------

.. autosummary::
    :toctree: generated/

    AtomGroup.ix
    AtomGroup.ix_array
    AtomGroup.n_atoms
    AtomGroup.n_residues
    AtomGroup.n_segments
    AtomGroup.residues
    AtomGroup.segments
    AtomGroup.ts
    AtomGroup.unique
    AtomGroup.universe

Methods
-------

.. autosummary::
    :toctree: generated/

    AtomGroup.select_atoms
    AtomGroup.write

Analysis properties
-------------------

.. autosummary::
    :toctree: generated/

    AtomGroup.bbox
    AtomGroup.bsphere
    AtomGroup.dimensions
    AtomGroup.forces
    AtomGroup.positions
    AtomGroup.velocities


Topology methods
----------------

.. autosummary::
    :toctree: generated/

    AtomGroup.accumulate
    AtomGroup.groupby
    AtomGroup.guess_bonds


Analysis methods
----------------

.. autosummary::
    :toctree: generated/

    AtomGroup.center
    AtomGroup.center_of_geometry
    AtomGroup.center_of_mass
    AtomGroup.centroid

Transformation
--------------

.. autosummary::
    :toctree: generated/

    AtomGroup.pack_into_box
    AtomGroup.rotate
    AtomGroup.rotateby
    AtomGroup.transform
    AtomGroup.translate
    AtomGroup.unwrap
    AtomGroup.wrap

Topology objects
----------------

.. autosummary::
    :toctree: generated/

    AtomGroup.angle
    AtomGroup.bond
    AtomGroup.dihedral
    AtomGroup.improper