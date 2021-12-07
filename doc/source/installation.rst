.. -*- coding: utf-8 -*-

====================
Installation
====================

The latest versions of **MDAnalysis** can be installed using `conda` or `pip`. 
Currently, the conda releases only support serial calculations.
If you plan to use the parallel OpenMP algorithms, you need to 
install MDAnalysis with pip and have a working OpenMP installation.

MDAnalysis has a separate :ref:`test suite <mdanalysistests>` **MDAnalysisTests** that is required to run the test cases and examples. 
The test files change less frequently, take up around 90 MB of space, 
and are not needed for daily use of MDAnalysis. However, they are often used in examples,
including many in this User Guide. If you are not interested in developing 
MDAnalysis or using the example files, you most likely don't need the tests. If you want to 
run examples in the User Guide, install the tests. 
The tests are distributed separately from the main package.

.. note::

    If you are installing on Windows, you must have
    Microsoft Visual C++ 14.0 installed. If your installation
    fails with the error message:

        error: Microsoft Visual C++ 14.0 is required. Get it with "Build Tools for Visual Studio": https://visualstudio.microsoft.com/downloads/
    
    Try installing Build Tools for Visual Studio from
    https://visualstudio.microsoft.com/downloads/ (scroll
    down to the Tools for Visual Studio section).


If you encounter any issues following these instructions, please
ask for help on the `user mailing list`_.

.. _`user mailing list`: https://groups.google.com/forum/#!forum/mdnalysis-discussion

conda
=====
To install the latest stable version of MDAnalysis via ``conda``, use the following command. This installs all dependencies needed for full analysis functionality (excluding external programs such as `HOLE`_):

.. code-block:: bash

    conda install -c conda-forge mdanalysis

To upgrade:

.. code-block:: bash

    conda update mdanalysis

To install the tests:

.. code-block:: bash

    conda install -c conda-forge MDAnalysisTests

If you intend to use MDAnalysis in JupyterLab, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    conda install -c conda-forge nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

pip
=====
The following command will install or upgrade the latest stable version of MDAnalysis via ``pip``, with core dependencies. This means that some packages required by specific analysis modules will not be installed.

.. code-block:: bash

    pip install --upgrade MDAnalysis

If you need to install a fully-featured MDAnalysis, add the ``analysis`` tag. As with ``conda``, this will not install external programs such as `HOLE`_.

.. code-block:: bash

    pip install --upgrade MDAnalysis[analysis]

To install/upgrade tests:

.. code-block:: bash

    pip install --upgrade MDAnalysisTests

If you intend to use MDAnalysis in JupyterLab, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    pip install nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager


Development versions
====================
To install development versions of MDAnalysis, you can compile it from source. In order to install from source, you will need ``numpy`` and ``cython``. See :ref:`create-virtual-environment` for instructions on how to create a full development environment.

.. code-block:: bash

    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    # assuming you have already installed required dependencies
    pip install -e package/

And to install the test suite:

.. code-block:: bash

    pip install -e testsuite/


Testing
-------

The tests rely on the `pytest` and `numpy` packages, which must also be installed. Run tests with: 

.. code-block:: bash

    pytest --disable-pytest-warnings --pyargs MDAnalysisTests

All tests should pass (i.e. no FAIL, ERROR); SKIPPED or XFAIL are ok. If anything fails or gives an error, 
`ask on the user mailing list <http://users.mdanalysis.org/>`_ or `raise an issue <https://github.com/MDAnalysis/mdanalysis/issues>`_.

Testing MDAnalysis can take a while, as there are quite a few tests. 
The plugin `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`_ can be used to run tests in parallel.

.. code-block:: bash

    pip install pytest-xdist
    pytest --disable-pytest-warnings --pyargs MDAnalysisTests --numprocesses 4

Custom compiler flags and optimised installations
-------------------------------------------------

You can pass any additional compiler flags for the C/C++ compiler using the ``extra_cflags`` variable in ``setup.cfg``.
This allows you to add any additional compiler options required for your architecture. 

This option can also be used to tune your MDAnalysis installation for your current architecture using `-march`, `-mtune`, `-mcpu` and related compiler flags.
The use of these compiler flags allows CPU architecture specific optimisations. An example for an x86_64 machine would be to change the

.. code-block:: YAML

    #extra_cflags =

line in `setup.cfg` to instead be:

.. code-block:: YAML

    extra_cflags = -march=native -mtune=native

Use of these flags can give a significant performance boost where the compiler can effectively autovectorise.

Be sure to use the recommended flags for your target architecture. For example, ARM platforms recommend using ``-mcpu`` *instead* of ``-mcpu``, while
PowerPC platforms prefer *both* ``-mcpu`` and ``-mtune``.

Full dicussion of the these flags is available elsewhere and a list of supported options should be provided by your compiler. The list for GCC_ is provided here. 

.. warning::
    Use of these compiler options is considered **advanced** and may reduce the binary compatibility of MDAnalysis significantly, especially if using `-march`,
    making it usable only on a matching CPU architecture to the one it is compiled on. We **strongly** recommend that you run the test suite on your intended platform
    before proceeding with analysis.

For example, if you plan on using MDAnalysis on a heterogenous system, such as a supercomputer, where the login node you compile on and the
compute node you run MDAnalysis on have different CPU architectures, you should avoid this option unless you know what you are doing, or you should compile directly on the compute node.

Additional datasets
===================

:ref:`MDAnalysisData is an additional package <mdanalysisdata>` with datasets that can be used in example tutorials. You can install it with ``conda`` or ``pip``:

.. code-block:: bash

    # conda
    conda install -c conda-forge mdanalysisdata
    # pip
    pip install --upgrade MDAnalysisData

This installation does not download all the datasets; instead, the datasets are cached when they are first downloaded using a Python command. 


.. _`HOLE`: http://www.holeprogram.org
.. _GCC: https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
.. _mdanalysisdata: https://www.mdanalysis.org/MDAnalysisData
