================================
  MDAnalysis User Guide
================================

|build|


MDAnalysis_ is a Python library for the analysis of computer simulations of many-body systems at the molecular scale, spanning use cases from interactions of drugs with proteins to novel materials. It is widely used in the scientific community and is written by scientists for scientists.

It works with a wide range of popular simulation packages including GROMACS, Amber, NAMD, CHARMM, DL_POLY, HOOMD, LAMMPS and many others — see the lists of supported `trajectory formats`_ and `topology formats`_.

This repository contains the `user guide`_ for MDAnalysis. The user guide originated as a `project for Google Season of Docs 2019-2020`_.

The code in the user guide is released under the `GNU LGPL v2.1+`_ license, while text and images are released under the `CC BY-SA 4.0 DEED`_ license.

.. _numfocus-fiscal-sponsor-attribution:

The MDAnalysis project uses an `open governance model`_ and is fiscally sponsored by `NumFOCUS`_. Consider making
a `tax-deductible donation`_ to help the project pay for developer time, professional services, travel, workshops, and a variety of other needs.

.. image:: https://www.mdanalysis.org/public/images/numfocus-sponsored-small.png
   :alt: NumFOCUS (Fiscally Sponsored Project)
   :target: https://numfocus.org/project/mdanalysis
   :align: center

This project is bound by a `Code of Conduct`_.

|powered_by_MDA|

If you use MDAnalysis_ in your project consider letting your users and the world know about it by displaying the MDAnalysis_ badge! `Embedding code`_ is available for different markups.

Contributing
============

If you would like to contribute documentation to the `user guide`_, please follow the instructions in the `contributing to the user guide`_ page. All contributions must comply with the MDAnalysis `AI policy`_.


Citation
========

When using MDAnalysis in published work, please cite the following
two papers:

*   R\. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy,
    M\. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski,
    S\. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis:
    A Python package for the rapid analysis of molecular
    dynamics simulations. In S. Benthall and S. Rostrup,
    editors, Proceedings of the 15th Python in Science
    Conference, pages 102-109, Austin, TX, 2016. SciPy.
    doi: `10.25080/Majora-629e541a-00e`_
*   N. Michaud-Agrawal, E. J. Denning, T. B. Woolf,
    and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular
    Dynamics Simulations. *J. Comput. Chem.* **32** (2011), 2319--2327.
    doi: `10.1002/jcc.21787`_

For citations of included algorithms and sub-modules please see the references_.

.. Footnotes

.. _`CC BY-SA 4.0 DEED`: https://creativecommons.org/licenses/by-sa/4.0/
.. _`GNU LGPL v2.1+`: https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html
.. _`user guide`: https://www.mdanalysis.org/UserGuide/index.html
.. _`many popular formats`: https://www.mdanalysis.org/UserGuide/formats/index.html
.. _MDAnalysis: https://www.mdanalysis.org
.. _LICENSE: https://github.com/MDAnalysis/mdanalysis/blob/develop/LICENSE
.. _`project for Google Season of Docs 2019-2020`: https://developers.google.com/season-of-docs/docs/2019/participants/project-mdanalysis
.. _NumFOCUS: https://numfocus.org/
.. _open governance model: https://www.mdanalysis.org/pages/governance/
.. _tax-deductible donation: https://numfocus.org/donate-to-mdanalysis
.. _`Code of Conduct`: https://www.mdanalysis.org/conduct/
.. _trajectory formats: https://docs.mdanalysis.org/documentation_pages/coordinates/init.html#id1
.. _topology formats: https://docs.mdanalysis.org/documentation_pages/topology/init.html#supported-topology-formats
.. _`AI policy`: https://github.com/MDAnalysis/mdanalysis/blob/develop/AI_POLICY.md
.. _`contributing to the user guide`: https://userguide.mdanalysis.org/stable/contributing_docs.html
.. _`10.1002/jcc.21787`: https://dx.doi.org/10.1002/jcc.21787
.. _`10.25080/Majora-629e541a-00e`: https://doi.org/10.25080/Majora-629e541a-00e
.. _references: https://docs.mdanalysis.org/documentation_pages/references.html
.. _Embedding code: https://www.mdanalysis.org/pages/citations/#powered-by-mdanalysis

.. |build| image:: https://github.com/MDAnalysis/UserGuide/actions/workflows/gh-ci.yaml/badge.svg
   :alt: Build Status
   :target: https://github.com/MDAnalysis/UserGuide/actions/workflows/gh-ci.yaml

.. |powered_by_MDA| image:: https://img.shields.io/badge/Powered%20by-MDAnalysis-orange.svg?logoWidth=15&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
   :alt: Powered by MDAnalysis
   :target: https://www.mdanalysis.org
