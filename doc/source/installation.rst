.. -*- coding: utf-8 -*-

====================
Installation Guide
====================

The latest version of **MDAnalysis** can be installed using :ref:`conda-installation` (recommended), 
or :ref:`pip-installation`. For development and contribution users can :ref:`source-installation`.
 
Currently, the `conda`_ releases only support serial calculations. If you plan to use 
the parallel `OpenMP`_ algorithms, you need to install MDAnalysis with `pip`_ and have a 
working `OpenMP`_ installation.

**MDAnalysisTests** is optional and it is a separate :ref:`test suite <mdanalysistests>` used 
for verifying MDAnalysis installations and running User Guide examples. It is not required for 
daily use, but some tutorials depend on it. The package is 90 MB and does not change frequently. 
If you plan to use tutorial examples or contribute to MDAnalysis development, installing it is recommended.

MDAnalysis requires Python 3.8 or later. Ensure you have an appropriate Python version installed before proceeding.

.. note::
    MDAnalysis supports **Linux**, **macOS**, and **Windows** (Python 3.8+). 

    - If you encounter errors on **Windows** related to **Microsoft Visual C++ 14.0**, install the required **Build Tools for Visual Studio** from: `Microsoft Visual Studio Downloads`_. 
    - If you encounter any other issues following these instructions, seek help on `GitHub Discussions (Installation)`_.    

.. warning::
    `conda`_ and `pip`_ installations **do not include** external programs such as `HOLE`_.

.. _conda-installation:

conda 
======
 
For most users, `conda`_ is the easiest way to install MDAnalysis because it manages dependencies effectively. 
We highly recommend creating a new environment for MDAnalysis to ensure a clean installation of MDAnalysis and 
its dependencies. We further recommend that you install and use `mamba`_, a faster drop-in replacement for `conda`_. 

If you don't have `conda`_ you can follow the `conda installation instructions`_. You can then create the environment 
and install `mamba`_ with: 

.. code-block:: bash

    conda create --name mdanalysis
    conda activate mdanalysis
    conda install -c conda-forge mamba

To install the latest stable version of MDAnalysis via `mamba`_ with all dependencies for full functionality, use the following command. 

.. code-block:: bash

    mamba install -c conda-forge mdanalysis

To upgrade use:

.. code-block:: bash

    mamba update mdanalysis

To install the :ref:`test suite <mdanalysistests>` use:

.. code-block:: bash

    mamba install -c conda-forge MDAnalysisTests

If you intend to use MDAnalysis in **JupyterLab**, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    conda install -c conda-forge nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

.. _pip-installation:

pip 
===

The following command will install or upgrade the latest stable version of MDAnalysis via `pip`_ with core dependencies. 
This means that some packages required by specific analysis modules will not be installed.

.. code-block:: bash

    pip install --upgrade MDAnalysis


If you need to install a fully-featured MDAnalysis, add the ``analysis`` tag. 

.. code-block:: bash

    pip install --upgrade MDAnalysis[analysis]

To install or upgrade the :ref:`test suite <mdanalysistests>`:

.. code-block:: bash

    pip install --upgrade MDAnalysisTests

If you intend to use MDAnalysis in **JupyterLab**, you will have to install
an extra package for the progress bar in analysis classes:

.. code-block:: bash

    pip install nodejs
    jupyter labextension install @jupyter-widgets/jupyterlab-manager

.. _source-installation:

Install from source
===================

If you plan to develop or contribute to MDAnalysis, follow these instructions first. Then, refer to :ref:`contributing` for additional guidelines.

To install the latest development versions of MDAnalysis (including unreleased features and bug fixes) you can compile it from `source`_. Use the following commands:
   
.. code-block:: bash

    # Clone the MDAnalysis repository
    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    # assuming you have already installed required dependencies
    pip install -e package/

.. note::
    To install from `source`_, you will need ``numpy`` and ``cython``. See :ref:`create-virtual-environment` for instructions on how to create a full development environment.

To install the test suite:

.. code-block:: bash
    
    # Install the test suite (optional but recommended for contributors)
    pip install -e testsuite/

Testing
-------

The tests rely on the `pytest`_ and `numpy`_ packages, which must also be installed. Run tests with:

.. code-block:: bash

    pytest --disable-pytest-warnings --pyargs MDAnalysisTests

All tests should pass (i.e. no FAIL, ERROR); SKIPPED or XFAIL are ok. If anything fails or gives an error,
ask on `GitHub Discussions`_ or `raise an issue`_.

Running MDAnalysis tests can take some time, as there are many test cases.
The plugin `pytest-xdist`_ can be used to run tests in parallel.

.. code-block:: bash

    pip install pytest-xdist
    pytest --disable-pytest-warnings --pyargs MDAnalysisTests --numprocesses 4

Custom compiler flags and optimised installations
-------------------------------------------------

You can pass any additional compiler flags for the C/C++ compiler using the ``extra_cflags`` variable in ``setup.cfg``.
This allows you to add any additional compiler options required for your architecture.

For example, ``extra_cflags`` can be used to tune your MDAnalysis installation for your current architecture using the `-march`, `-mtune`, `-mcpu` and related compiler flags.
The optimal compiler flags depend on your CPU architecture. An example for an x86_64 machine would be to change the line in `setup.cfg` as follows:

.. code-block:: diff

	- #extra_cflags =
	+ extra_cflags = -march=native -mtune=native

Use of these flags can give a significant performance boost where the compiler can effectively autovectorise.

Be sure to use the recommended flags for your target architecture. For example, ARM platforms recommend using ``-mcpu`` instead of ``-march``, while 
PowerPC platforms prefer both ``-mcpu`` and ``-mtune``.

Full discussion of these flags is available elsewhere (such as here in this `wiki`_ or in this `ARM`_ blog post) and a list of supported options should be provided by your compiler. The list for GCC_ is provided here.

.. warning::
    These compiler options are **advanced** and may reduce binary compatibility. In particular, using `-march` may restrict MDAnalysis to the exact CPU architecture it was compiled on.
    We **strongly recommend** that you run the test suite on your intended platform before proceeding with analysis.

In cases where you might encounter multiple CPU architectures (e.g. on a supercomputer where the login node and compute node have different architectures), you should avoid changing these options unless you are experienced with compiling software in these situations.

Additional datasets
===================

MDAnalysisData_ is an additional package with datasets that can be used in example tutorials. You can install it with `conda`_ or `pip`_:

.. code-block:: bash

    # conda
    conda install -c conda-forge mdanalysisdata
    # pip
    pip install --upgrade MDAnalysisData

This installation does not download all the datasets; instead, the datasets are cached when they are first downloaded using a Python command.


.. _HOLE: http://www.holeprogram.org 
.. _GCC: https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
.. _MDAnalysisData: https://www.mdanalysis.org/MDAnalysisData/
.. _wiki: https://wiki.gentoo.org/wiki/GCC_optimization#-march
.. _ARM: https://community.arm.com/arm-community-blogs/b/tools-software-ides-blog/posts/compiler-flags-across-architectures-march-mtune-and-mcpu
.. _pip: https://pip.pypa.io/en/latest/index.html 
.. _mamba: https://anaconda.org/conda-forge/mamba
.. _conda: https://docs.conda.io/projects/conda/en/latest/
.. _source: https://github.com/MDAnalysis/mdanalysis/
.. _GitHub Discussions (Installation): https://github.com/MDAnalysis/mdanalysis/discussions/categories/installation
.. _Microsoft Visual Studio Downloads: https://visualstudio.microsoft.com/downloads/
.. _pytest: https://docs.pytest.org/en/stable/
.. _numpy: https://numpy.org/
.. _Github discussions: https://github.com/MDAnalysis/mdanalysis/discussions 
.. _raise an issue: https://github.com/MDAnalysis/mdanalysis/issues
.. _pytest-xdist: https://github.com/pytest-dev/pytest-xdist 
.. _OpenMP: https://www.openmp.org/
.. _conda installation instructions: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html