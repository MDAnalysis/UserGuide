.. -*- coding: utf-8 -*-

====================
Installation
====================

The latest versions of MDAnalysis can be installed using `conda` or `pip`. 
Currently, the conda releases only support serial calculations.
If you plan to use the parallel OpenMP algorithms, you need to 
install MDAnalysis with pip and have a working OpenMP installation.

MDAnalysis has a separate test suite that is required to run the test cases and examples. 
The test files change less frequently, take up around 90 MB of space, 
and are not needed for daily use of MDAnalysis.
They are distributed separately from the main package. 

conda
=====
To install the latest stable version of MDAnalysis via Anaconda, use the following command:

.. code-block:: bash

    conda install -c conda-forge mdanalysis

To upgrade:

.. code-block:: bash

    conda update mdanalysis

To install the tests:

.. code-block:: bash

    conda install -c conda-forge MDAnalysisTests

pip
=====
To install the latest stable version of MDAnalysis via pip, use the following command:

.. code-block:: bash

    pip install --upgrade MDAnalysis

To upgrade:

.. code-block:: bash

    conda update mdanalysis

To install tests:

.. code-block:: bash

    pip install --upgrade MDAnalysisTests


Development versions
====================
To install development versions of MDAnalysis, you can compile it from source.

.. code-block:: bash

    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    pip install -e .
