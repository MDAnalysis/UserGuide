.. -*- coding: utf-8 -*-
.. _standard-selections:

==========================================
Standard residues in MDAnalysis selections
==========================================

.. _protein-selection:

Proteins
========
The residue names listed here are accessible via the "protein" keyword in the :ref:`selections`. 

The below names are drawn from the CHARMM 27, OPLS-AA, GROMOS 53A6, AMBER 03, and AMBER 99sb*-ILDN force fields.

+------+------+------+------+------+------+------+------+
| ACE  | CASF | CLYS | CYS2 | HIS1 | LYN  | NGLY | NTRP |
+------+------+------+------+------+------+------+------+
| ALA  | CASN | CME  | CYSH | HIS2 | LYS  | NHID | NTYR |
+------+------+------+------+------+------+------+------+
| ALAD | CASP | CMET | CYX  | HISA | LYSH | NHIE | NVAL |
+------+------+------+------+------+------+------+------+
| ARG  | CCYS | CPHE | DAB  | HISB | MET  | NHIP | ORN  |
+------+------+------+------+------+------+------+------+
| ARGN | CCYX | CPRO | GLH  | HISD | MSE  | NILE | PGLU |
+------+------+------+------+------+------+------+------+
| ASF  | CGLN | CSER | GLN  | HISE | NALA | NLEU | PHE  |
+------+------+------+------+------+------+------+------+
| ASH  | CGLU | CTHR | GLU  | HISH | NARG | NLYS | PRO  |
+------+------+------+------+------+------+------+------+
| ASN  | CGLY | CTRP | GLUH | HSD  | NASN | NME  | QLN  |
+------+------+------+------+------+------+------+------+
| ASN1 | CHID | CTYR | GLY  | HSE  | NASP | NMET | SER  |
+------+------+------+------+------+------+------+------+
| ASP  | CHIE | CVAL | HID  | HSP  | NCYS | NPHE | THR  |
+------+------+------+------+------+------+------+------+
| ASPH | CHIP | CYM  | HIE  | HYP  | NCYX | NPRO | TRP  |
+------+------+------+------+------+------+------+------+
| CALA | CILE | CYS  | HIP  | ILE  | NGLN | NSER | TYR  |
+------+------+------+------+------+------+------+------+
| CARG | CLEU | CYS1 | HIS  | LEU  | NGLU | NTHR | VAL  |
+------+------+------+------+------+------+------+------+

----------------
Protein backbone
----------------

Protein backbone atoms in MDAnalysis belong to a recognised protein residue and have the atom names CA, C, O, N. 

.. _nucleic-selection:

Nucleic acids
=============

The residue names listed here are accessible via the "nucleic" keyword in the :ref:`selections`. 

The below names are drawn from largely from the CHARMM force field.

+-----+----+----+---+-----+-----+-----+-----+
| ADE | DA | RU | T | DA5 | DA3 | RA5 | RA3 |
+-----+----+----+---+-----+-----+-----+-----+
| URA | DC | RG | U | DC5 | DC3 | RU5 | RU3 |
+-----+----+----+---+-----+-----+-----+-----+
| CYT | DG | RC | C | DG5 | DG3 | RG5 | RG3 |
+-----+----+----+---+-----+-----+-----+-----+
| GUA | DT | A  | G | DT5 | DT3 | RC5 | RC3 |
+-----+----+----+---+-----+-----+-----+-----+
| THY | RA |    |   |     |     |     |     |
+-----+----+----+---+-----+-----+-----+-----+

----------------
Nucleic backbone
----------------

Nucleic backbone atoms in MDAnalysis belong to a recognised nucleic acid residue and have the atom names P, C5', C3', O3', O5'.

.. _nucleobase-selection:

-----------
Nucleobases
-----------

Nucleobase atoms from nucleic acid residues are recognised based on their names in CHARMM.

+----+----+----+----+----+----+----+-----+
| N9 | C8 | C4 | C2 | C6 | N2 | O2 | O4  |
+----+----+----+----+----+----+----+-----+
| N7 | C5 | N3 | N1 | O6 | N6 | N4 | C5M |
+----+----+----+----+----+----+----+-----+

--------------
Nucleic sugars
--------------

Nucleic sugar atoms from nucleic acid residues are recognised by MDAnalysis if they have the atom names C1', C2', C3', C4', O2', O4', O3'. 