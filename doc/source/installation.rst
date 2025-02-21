.. -*- coding: utf-8 -*-

====================
Installation Guide
====================

**MDAnalysis** can be installed using `conda <https://docs.conda.io/projects/conda/en/latest/>`_ (recommended), `pip <https://pip.pypa.io/en/latest/index.html>`_, or from `source <https://github.com/MDAnalysis/mdanalysis/>`_ for development.

.. note::
    MDAnalysis supports **Linux**, **macOS**, and **Windows** (Python 3.8+).
    If you encounter installation issues, seek help on `GitHub Discussions (Installation) <https://github.com/MDAnalysis/mdanalysis/discussions/categories/installation>`_.

Conda installation
=====

For most users, `conda <https://docs.conda.io/projects/conda/en/latest/>`_ is the easiest way to install MDAnalysis. We highly recommend creating a new environment for MDAnalysis.
This will ensure that you have a clean installation of MDAnalysis and its dependencies, and will not interfere with other packages you may have installed.

.. code-block:: bash

    conda create --name mdanalysis
    conda activate mdanalysis
    conda install -c conda-forge mdanalysis

To upgrade to the latest version:

.. code-block:: bash

    conda update mdanalysis

To install the **test suite** (useful for examples and validation):

.. code-block:: bash

    conda install MDAnalysisTests

.. warning::
    Conda installations **do not support OpenMP**. 
    If you need **parallel OpenMP calculations**, install MDAnalysis using `pip <https://pip.pypa.io/en/latest/index.html>`_.


If you intend to use MDAnalysis in JupyterLab, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    conda install -c conda-forge nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

Pip installation
=====

If you do not use Conda or require **OpenMP support**, install MDAnalysis via `pip <https://pip.pypa.io/en/latest/index.html>`_:

.. code-block:: bash

    pip install --upgrade MDAnalysis

For full functionality with **analysis modules**:

.. code-block:: bash

    pip install --upgrade MDAnalysis[analysis]

To install the **test suite**:

.. code-block:: bash

    pip install --upgrade MDAnalysisTests

If you intend to use MDAnalysis in JupyterLab, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    pip install nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

.. note::
    If you are installing on **Windows**, and encounter errors related to **Microsoft Visual C++ 14.0**, install **Build Tools for Visual Studio** from `here <https://visualstudio.microsoft.com/downloads/>`_.

----------------------

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
`ask on GitHub Discussions <https://github.com/MDAnalysis/mdanalysis/discussions>`_ or `raise an issue <https://github.com/MDAnalysis/mdanalysis/issues>`_.

Testing MDAnalysis can take a while, as there are quite a few tests.
The plugin `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`_ can be used to run tests in parallel.

.. code-block:: bash

    pip install pytest-xdist
    pytest --disable-pytest-warnings --pyargs MDAnalysisTests --numprocesses 4

Custom compiler flags and optimised installations
-------------------------------------------------

You can pass any additional compiler flags for the C/C++ compiler using the ``extra_cflags`` variable in ``setup.cfg``.
This allows you to add any additional compiler options required for your architecture.

For example, ``extra_cflags`` can be used to tune your MDAnalysis installation for your current architecture using the `-march`, `-mtune`, `-mcpu` and related compiler flags.
*Which* particular compiler flags to use depends on your CPU architecture. An example for an x86_64 machine would be to change the line in `setup.cfg` as follows:

.. code-block:: diff

	- #extra_cflags =
	+ extra_cflags = -march=native -mtune=native

Use of these flags can give a significant performance boost where the compiler can effectively autovectorise.

Be sure to use the recommended flags for your target architecture. For example, ARM platforms recommend using ``-mcpu`` *instead* of ``-mcpu``, while
PowerPC platforms prefer *both* ``-mcpu`` and ``-mtune``.

Full dicussion of the these flags is available elsewhere (such as here in this wiki_ or in this ARM_ blog post) and a list of supported options should be provided by your compiler. The list for GCC_ is provided here.

.. warning::
    Use of these compiler options is considered **advanced** and may reduce the binary compatibility of MDAnalysis significantly, especially if using `-march`,
    making it usable only on a matching CPU architecture to the one it is compiled on. We **strongly** recommend that you run the test suite on your intended platform
    before proceeding with analysis.

In cases where you might encounter multiple CPU architectures (e.g. on a supercomputer where the login node and compute node have different architectures), you should avoid changing these options unless you are experienced with compiling software in these situations.

Additional datasets
===================

MDAnalysisData_ is an additional package with datasets that can be used in example tutorials. You can install it with ``conda`` or ``pip``:

.. code-block:: bash

    # conda
    conda install -c conda-forge mdanalysisdata
    # pip
    pip install --upgrade MDAnalysisData

This installation does not download all the datasets; instead, the datasets are cached when they are first downloaded using a Python command.


.. _`HOLE`: http://www.holeprogram.org
.. _GCC: https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
.. _MDAnalysisData: https://www.mdanalysis.org/MDAnalysisData/
.. _wiki: https://wiki.gentoo.org/wiki/GCC_optimization#-march
.. _ARM: https://community.arm.com/arm-community-blogs/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
