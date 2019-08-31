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
and are not needed for daily use of MDAnalysis. If you are not interested in developing 
MDAnalysis or using the example files, you most likely don't need the tests. 
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

Testing
-------

The tests rely on the `pytest` and `numpy` packages, which must also be installed. Run tests with: 

.. code-block:: bash

    pytest --disable-pytest-warnings --pyargs MDAnalysisTests

All tests should pass (i.e. no FAIL, ERROR); SKIPPED or XFAIL are ok. If anything fails or gives an error, 
[ask on the user mailing list](http://users.mdanalysis.org/) or [raise an issue](https://github.com/MDAnalysis/mdanalysis/issues).

Testing MDAnalysis can take a while, as there are quite a few tests. 
The plugin [pytest-xdist](https://github.com/pytest-dev/pytest-xdist) can be used to run tests in parallel.

.. code-block:: bash

    pip install pytest-xdist
    pytest --disable-pytest-warnings --pyargs MDAnalysisTests --numprocesses 4
