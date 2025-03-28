.. -*- coding: utf-8 -*-

====================
Installation Guide
====================

The latest version of **MDAnalysis** can be installed using :ref:`mamba-installation` (recommended),
or :ref:`pip-installation`. If you need to install MDAnalysis from source, follow :ref:`source-installation`
for stable releases. If you are a contributor or developer working on MDAnalysis, follow the
:ref:`development-installation` for editable installation and code modifications.

Currently, the :ref:`mamba-installation` installs a version of MDAnalysis that does **not** include `OpenMP`_
acceleration due to limitations in precompiled conda packages. If you need `OpenMP`_-enabled features,
consider installing via :ref:`pip-installation` and compiling from source with a compatible compiler and
OpenMP support.

**MDAnalysisTests** is optional and is a separate :ref:`test suite <mdanalysistests>` used
for verifying MDAnalysis installations and running User Guide examples. It is not required for
daily use, but some tutorials rely on it. The package is approximately 90 MB and is not updated frequently.
If you plan to use tutorial examples or contribute to MDAnalysis, installing it is recommended.

MDAnalysis supports specific Python versions. For the latest supported versions, please refer to the
`CHANGELOG`_ or the `PyPI`_. Ensure you have a compatible Python version installed before proceeding.

.. note::
    MDAnalysis supports **Linux**, **macOS**, and **Windows**.

    - If you encounter errors on **Windows** related to **Microsoft Visual C++ 14.0**, install the required **Build Tools for Visual Studio** from: `Microsoft Visual Studio Downloads`_.
    - If you encounter any other issues following these instructions, seek help on `GitHub Discussions (Installation)`_.

.. _mamba-installation:

mamba/conda
===========

For most users, `mamba`_ is the **recommended** way to install MDAnalysis. It is a faster drop-in
replacement for `conda`_ and efficiently handles dependencies. We also recommend creating a new
environment for MDAnalysis to ensure a clean and isolated installation.

MDAnalysis supports two common ways to install with `mamba`_:

**1. Recommended: Use Mambaforge**

`Mambaforge`_ is a minimal `conda`_ distribution that includes `mamba`_ by default, uses the `conda-forge`_ channel
(preferred for MDAnalysis) and avoids clutter from `Anaconda`_ base packages.

To install MDAnalysis with `Mambaforge`_, follow the `mambaforge installation instructions`_ and then
run the following commands to install MDAnalysis:

.. code-block:: bash

    mamba create -n mdanalysis -c conda-forge mdanalysis
    mamba activate mdanalysis

To install the :ref:`test suite <mdanalysistests>` (optional, ~90 MB):

.. code-block:: bash

    mamba install -c conda-forge MDAnalysisTests

To upgrade, use:

.. code-block:: bash

    mamba update mdanalysis

If you're using MDAnalysis in **JupyterLab**, install `ipywidgets`_ for progress bar support:

.. code-block:: bash

    mamba install -c conda-forge ipywidgets

**2. Alternative: Use conda with conda-forge**

If you prefer using `conda`_ directly (e.g., installed via `Anaconda installer`_), create a new environment with packages from the `conda-forge`_ channel:

.. code-block:: bash

    conda create -n mdanalysis -c conda-forge mdanalysis
    conda activate mdanalysis

To install optional packages:

.. code-block:: bash

    conda install -c conda-forge MDAnalysisTests ipywidgets

.. _pip-installation:

pip
===

The following command will install or upgrade the latest stable version of MDAnalysis via `pip`_ with core dependencies.
This means that some packages required by specific analysis modules will not be installed.

.. code-block:: bash

    pip install --upgrade MDAnalysis

MDAnalysis offers several optional dependency groups (“extra tags”) that you can install with pip to enable additional features:

- ``analysis``: optional dependencies for various analysis modules.
- ``extra_formats``: support for additional file formats.
- ``parallel``: install dask to enable parallel analysis across multiple cores or machines.
- ``doc``: install documentation build dependencies (e.g. Sphinx, themes).

To see the full list of available extras, check the `pyproject.toml file`_.

To install with optional features, use pip with one or more extras:

.. code-block:: bash

    pip install --upgrade MDAnalysis[analysis,extra_formats,parallel]

To install or upgrade the :ref:`test suite <mdanalysistests>` (optional, ~90 MB):

.. code-block:: bash

    pip install --upgrade MDAnalysisTests

If you're using MDAnalysis in **JupyterLab**, install `ipywidgets`_ for progress bar support:

.. code-block:: bash

    pip install ipywidgets

.. _source-installation:

Install from source
===================

This section is for **users** who want to install a **specific stable version** of MDAnalysis
(e.g., for reproducibility or to debug a release). This is **not** intended for development or code contributions.
If you're a developer, see :ref:`development-installation`.

To install from source, you must first ensure that your environment already contains all necessary dependencies.
We recommend you check the :ref:`create-virtual-environment` to set up a clean development environment
with all required dependencies. We recommend using `mamba`_ or `conda`_ to manage this setup.

Follow the steps in:

- `Set up with conda-forge`_
- `Set up with mamba`_

MDAnalysis uses different Git branches for different purposes:

- **package-X.Y.Z**: These branches contain specific stable releases (e.g., ``package-2.9.0``) and are intended for users who need reproducibility, debugging, or version locking. You can browse them on the GitHub `Branches page`_ or view official releases on the `Releases page`_.
- **develop**: The `develop`_ branch is the default and active development branch. If you're contributing to MDAnalysis or want the latest features before they're released, base your work on this branch (see :ref:`development-installation`).

To install the **latest stable release** from source:

.. code-block:: bash

    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    git checkout release-<latest-version>
    pip install .

To run tests:

.. code-block:: bash

    pip install MDAnalysisTests

.. _development-installation:

Install for development
==========================

This section is for **contributors and developers** who want to modify MDAnalysis
or contribute to its development. You’ll install MDAnalysis in **editable mode**, which allows you to make code changes
without reinstalling the package each time.

Before installing, follow the steps described in the :ref:`source-installation` to set up a clean environment
with all dependencies.

Once your environment is ready, install MDAnalysis in editable mode:

.. code-block:: bash

    git clone https://github.com/MDAnalysis/mdanalysis
    cd mdanalysis
    git switch develop    # work on the develop branch
    git checkout -b my-feature-branch  # Create a new branch for your changes
    pip install -e package/

Installing in editable mode (`-e`) means that any changes you make to the MDAnalysis source code are immediately
available without needing to reinstall.

If you plan to run tests you can install the test suite:

.. code-block:: bash

    # Install the test suite (optional but recommended for contributors)
    pip install -e testsuite/

For more on how to set up and contribute, refer to the :ref:`contributing` guide.

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
The optimal compiler flags depend on your CPU architecture. Commonly recommended settings are:

.. code-block:: diff

	- #extra_cflags =
	+ extra_cflags = -march=native -mtune=native

Use of these flags can give a significant performance boost where the compiler can effectively autovectorise.

Be sure to use the recommended flags for your target architecture. In many cases,
``-march`` and ``-mcpu`` can be used together or separately, and their behavior depends
on your platform and compiler version. For example, ARM platforms often prefer ``-mcpu``,
while PowerPC may use both ``-mcpu`` and ``-mtune``. Refer to the `GCC documentation on -march`_
and `-mcpu`_ for platform-specific guidance.

Full discussion of these flags is available elsewhere (such as here in this `wiki`_ or in this `ARM`_ blog post) and a list of supported options should be provided by your compiler. The list for GCC_ is provided here.

.. warning::
    These compiler options are **advanced** and may reduce binary compatibility. In particular, using `-march` may restrict MDAnalysis to the exact CPU architecture it was compiled on.
    We **strongly recommend** that you run the test suite on your intended platform before proceeding with analysis.

In cases where you might encounter multiple CPU architectures (e.g. on a supercomputer where the login node and compute node have different architectures), you should avoid changing these options unless you are experienced with compiling software in these situations.

Additional datasets
===================

MDAnalysisData_ is an additional package with datasets that can be used in example tutorials. You can install it with `mamba`_ / `conda`_ or `pip`_:

.. code-block:: bash

    # conda
    conda install -c conda-forge mdanalysisdata
    # pip
    pip install --upgrade MDAnalysisData

This installation does not download all the datasets; instead, the datasets are cached when they are first downloaded using a Python command.

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
.. _mambaforge installation instructions: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _conda-forge: https://conda-forge.org/
.. _Anaconda: https://www.anaconda.com/products/distribution
.. _ipywidgets: https://ipywidgets.readthedocs.io/en/stable/
.. _Releases page: https://github.com/MDAnalysis/mdanalysis/releases
.. _Set up with conda-forge: https://userguide.mdanalysis.org/stable/contributing_code.html#with-conda-forge-packages
.. _Set up with mamba: https://userguide.mdanalysis.org/stable/contributing_code.html#with-mamba-conda
.. _develop: https://github.com/MDAnalysis/mdanalysis/tree/develop
.. _Branches page: https://github.com/MDAnalysis/mdanalysis/branches/all
.. _CHANGELOG: https://github.com/MDAnalysis/mdanalysis/blob/develop/package/CHANGELOG
.. _PyPI: https://pypi.org/project/MDAnalysis/
.. _pyproject.toml file: https://github.com/MDAnalysis/mdanalysis/blob/35d9d2e3ab08e7e6741b57fe02a7215fe3b91a6c/package/pyproject.toml#L69
.. _GCC documentation on -march: https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
.. _-mcpu: https://gcc.gnu.org/onlinedocs/gcc/ARM-Options.html
.. _Anaconda installer: https://docs.anaconda.com/free/anaconda/install/