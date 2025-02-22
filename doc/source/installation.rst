.. -*- coding: utf-8 -*-

====================
Installation Guide
====================

The latest version of **MDAnalysis** can be installed using:
`pip <https://pip.pypa.io/en/latest/index.html>`_ (most common, for general users),
`conda <https://docs.conda.io/projects/conda/en/latest/>`_ (recommended for better dependency management),
or from `source <https://github.com/MDAnalysis/mdanalysis/>`_ (for development and contributing). Currently, 
the `conda <https://docs.conda.io/projects/conda/en/latest/>`_ releases only support serial calculations. If you plan to use the parallel OpenMP algorithms, 
you need to install MDAnalysis with `pip <https://pip.pypa.io/en/latest/index.html>`_ and have a working OpenMP installation.

**MDAnalysisTests** is a optional and it is a separate :ref:`test suite <mdanalysistests>` used for verifying MDAnalysis installations and running User Guide examples. 
It is not required for daily use, but some tutorials depend on it. The package is 90 MB and does not change frequently. If you plan to use tutorial examples 
or contribute to MDAnalysis development, installing it is recommended.

.. note::
    MDAnalysis supports **Linux**, **macOS**, and **Windows** (Python 3.8+). 

    - If you encounter errors on **Windows** related to **Microsoft Visual C++ 14.0**, install the required **Build Tools for Visual Studio** from: `Microsoft Visual Studio Downloads <https://visualstudio.microsoft.com/downloads/>`_. 
    - If you encounter any other issues following these instructions, seek help on `GitHub Discussions (Installation) <https://github.com/MDAnalysis/mdanalysis/discussions/categories/installation>`_.    

.. warning::
    `Conda <https://docs.conda.io/projects/conda/en/latest/>`_ and `pip <https://pip.pypa.io/en/latest/index.html>`_ installations **do not include** external programs such as `HOLE <https://www.holeprogram.org//>`_.

Conda installation
=====
 
For most users, `conda <https://docs.conda.io/projects/conda/en/latest/>`_ is the easiest way to install MDAnalysis because it manages dependencies effectively. 
We highly recommend creating a new environment for MDAnalysis. We further recommend that you install and use ``mamba``, a faster drop-in replacement for ``conda``.

.. code-block:: bash

    conda create --name mdanalysis
    conda activate mdanalysis
    conda install -c conda-forge mamba

To install the latest stable version of MDAnalysis via ``conda``, use the following command. This installs all dependencies needed for full analysis functionality.

.. code-block:: bash

    mamba install -c conda-forge mdanalysis

To upgrade use:

.. code-block:: bash

    mamba update mdanalysis

To install the `test suite <https://userguide.mdanalysis.org/stable/datasets.html#mdanalysistests>`_ use:

.. code-block:: bash

    mamba install -c conda-forge MDAnalysisTests

If you intend to use MDAnalysis in **JupyterLab**, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    conda install -c conda-forge nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

.. warning::
    MDAnalysis installed via conda does not support OpenMP. If you need **parallel OpenMP calculations**, 
    **install MDAnalysis** using `pip <https://pip.pypa.io/en/latest/index.html>`_.

Pip installation
=====

The following command will install or upgrade the latest stable version of MDAnalysis via `pip <https://pip.pypa.io/en/latest/index.html>`_ with core dependencies. This means that some packages required by specific analysis modules will not be installed.

.. code-block:: bash

    pip install --upgrade MDAnalysis


If you need to install a fully-featured MDAnalysis, add the ``analysis`` tag. 

.. code-block:: bash

    pip install --upgrade MDAnalysis[analysis]

To install or upgrade the `test suite <https://userguide.mdanalysis.org/stable/datasets.html#mdanalysistests>`_:

.. code-block:: bash

    pip install --upgrade MDAnalysisTests

If you intend to use MDAnalysis in **JupyterLab**, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    pip install nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

Install from `source <https://github.com/MDAnalysis/mdanalysis/>`_
====================

To install the latest development versions of MDAnalysis (including unreleased features and bug fixes) you can compile it from `source <https://github.com/MDAnalysis/mdanalysis/>`_. Use the following commands:
   
.. code-block:: bash

    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    # assuming you have already installed required dependencies
    pip install -e package/

.. note::
    To install from `source <https://github.com/MDAnalysis/mdanalysis/>`_, you will need ``numpy`` and ``cython``. See :ref:`create-virtual-environment` for instructions on how to create a full development environment.

To install the test suite:

.. code-block:: bash

    pip install -e testsuite/

Testing
-------

The tests rely on the `pytest <https://docs.pytest.org/en/stable/>`_ and `numpy <https://numpy.org/>`_ packages, which must also be installed. Run tests with:

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
    These compiler options are **advanced** and may reduce binary compatibility. In particular, using `-march` may restrict MDAnalysis to the exact CPU architecture it was compiled on.
    We **strongly recommend** that you run the test suite on your intended platform before proceeding with analysis.

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
