.. -*- coding: utf-8 -*-
.. _PDB-label:

=======================================
PDB, ENT (Standard PDB file)
=======================================

.. include:: classes/PDB.txt

Reading in
==========

PDBReader that reads a `PDB-formatted`_ file, no frills.

The following *PDB records* are parsed (see `PDB coordinate section`_ for
details):

    - *CRYST1* for unitcell A,B,C, alpha,beta,gamma
    - *ATOM* or *HETATM* for serial,name,resName,chainID,resSeq,x,y,z,occupancy,tempFactor
    - *HEADER* (:attr:`header`), *TITLE* (:attr:`title`), *COMPND*
    (:attr:`compound`), *REMARK* (:attr:`remarks`)
    - all other lines are ignored

Reads multi-`MODEL`_ PDB files as trajectories.

.. _PDB-formatted:
    http://www.wwpdb.org/documentation/file-format-content/format32/v3.2.html
.. _PDB coordinate section:
    http://www.wwpdb.org/documentation/file-format-content/format32/sect9.html
.. _MODEL:
    http://www.wwpdb.org/documentation/file-format-content/format32/sect9.html#MODEL

=============  ============  ===========  =============================================
COLUMNS        DATA  TYPE    FIELD        DEFINITION
=============  ============  ===========  =============================================
1 -  6         Record name   "CRYST1"
7 - 15         Real(9.3)     a              a (Angstroms).
16 - 24        Real(9.3)     b              b (Angstroms).
25 - 33        Real(9.3)     c              c (Angstroms).
34 - 40        Real(7.2)     alpha          alpha (degrees).
41 - 47        Real(7.2)     beta           beta (degrees).
48 - 54        Real(7.2)     gamma          gamma (degrees).

1 -  6         Record name   "ATOM  "
7 - 11         Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 21        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
67 - 76        String        segID        (unofficial CHARMM extension ?)
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.
=============  ============  ===========  =============================================
