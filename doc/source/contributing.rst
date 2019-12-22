.. -*- coding: utf-8 -*-
.. _contributing:

===========================
Contributing to MDAnalysis
===========================

MDAnalysis is a free and open source project. It evolves and grows with the demands of its user base. 
The development team very much welcomes contributions from its users. 
Contributions can take many forms, such as:

    * **bug reports** or **enhancement requests** filed through the `Issue Tracker`_
    * **source code** for fixing bugs or new functionality (e.g. new analysis tasks)
    * **documentation** (either in the code or on the user guide)
    * **questions** or **discussions** on the `mdnalysis-discussion`_ mailing list

The MDAnalysis community subscribes to a `Code of Conduct`_ that all community
members agree and adhere to --- please read it.




Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

If you are brand new to MDAnalysis or open-source development, we recommend going 
through this page. If you are new to Git and version control, have a look at 
:ref:`version-control`. Get started working with the code of 
:ref:`the main MDAnalysis codebase <working-with-mdanalysis-code>` or 
:ref:`the user guide <


.. _version-control:

Version control, Git, and GitHub
================================

`Git <http://git-scm.com/>`_ is a version control system that allows many people to work together 
on a project. 
Working with Git can be one of the more daunting aspects of contributing 
to MDAnalysis.  Sticking to the guidelines below will help keep the process 
straightforward and mostly trouble free.  As always,
if you are having difficulties please feel free to ask for help.

The code is hosted on `GitHub <https://www.github.com/pydata/xarray>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. In short, the development workflow looks like this:

    #. :ref:`Fork the MDAnalysis repository from the mdanalysis account into your own account <forking-code-repo>`
    #. :ref:`Set up an isolated virtual environment for code development <create-virtual-environment>`
    #. :ref:`Build development versions of MDAnalysis and MDAnalysisTests on your computer into the virtual environment <build-mdanalysis-develop>`
    #. :ref:`Create a new branch off the develop branch <create-code-branch>`
    #. Add your new feature or bug fix
    #. Add and run tests 

Some great resources for learning Git:

    * the `GitHub help pages <http://help.github.com/>`_.
    * the `NumPy's documentation <http://docs.scipy.org/doc/numpy/dev/index.html>`_.
    * Matthew Brett's `Pydagogue <http://matthew-brett.github.com/pydagogue/>`_.

------------------------
Getting started with Git
------------------------

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _working-with-mdanalysis-code:

Working with the code
=====================

.. _forking-code-repo:

-------
Forking
-------

You will need your own fork to work on the code. Go to the `MDAnalysis project page <https://github.com/MDAnalysis/mdanalysis>`_ and hit the :guilabel:`Fork` button. You will
want to clone your fork to your machine:

.. code-block:: bash

    git clone https://github.com/your-user-name/mdanalysis.git
    cd mdanalysis
    git remote add upstream https://github.com/MDAnalysis/mdanalysis

This creates the directory `mdanalysis` and connects your repository to
the upstream (main project) MDAnalysis repository.

.. _create-virtual-environment:

----------------------------------
Creating a development environment
----------------------------------

To change code and test changes, you'll need to build both **MDAnalysis** and **MDAnalysisTests** 
from source. This requires a Python environment. We highly recommend that you use 
virtual environments. This allows you to have multiple experimental development versions 
of MDAnalysis which do not interfere with each other or your own stable version. 
Since MDAnalysis is split into the actual package and a test suite, you need to install 
both modules in development mode.

You can do this either with :ref:`conda <dev-with-conda>` or :ref:`pip <dev-with-pip>`.

.. _dev-with-conda:

With conda
----------

Install either `Anaconda <https://www.anaconda.com/download/>`_ 
or `miniconda <https://conda.io/miniconda.html>`_.
Make sure your conda is up to date:

.. code-block:: bash

    conda update conda

Create a new environment with ``conda create``. This will allow you to change code in 
an isolated environment without touching your base Python installation, and without 
touching existing environments that may have stable versions of MDAnalysis. :

.. code-block:: bash

    conda create --name mdanalysis-dev

Activate the environment to build MDAnalysis into it:

.. code-block:: bash

    conda activate mdanalysis-dev

To view your environments:

.. code-block:: bash

    conda info -e

To list the packages installed in your current environment:

.. code-block:: bash

    conda list

To return to your root environment:

.. code-block:: bash

    conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

.. _dev-with-pip:

With pip and virtualenv
-----------------------

Like conda, virtual environments managed with virtualenv allow you to use different versions of python and python packages for your different project. Unlike conda, virtualenv is not a general-purpose package manager; it leverages what is available on your system, and let you install python packages using pip.

To use virtual environments you have to install the virtualenv package first. This can be done with either pip or the package manager of your system:

.. code-block:: bash

    pip install virtualenv
    # or on ubuntu
    sudo apt install virtualenv
    # or on fedora
    sudo dnf install python-virtualenv

Virtual environments can be created per project directory.

.. code-block:: bash

    cd my-project/
    virtualenv my-project-env

This will create a new folder ``my-project-env``. This folder contains the virtual environment and all packages you have installed in it. To activate it in the current terminal run::

    source myproject-env/bin/activate

Now you can install packages via pip without affecting your global environment. The packages you install when the environment is activated will be available in terminal sessions that have the environment activated. You can deactivate the virtual environment by running::

    deactivate

The ``virtualenvwrapper`` `package <https://virtualenvwrapper.readthedocs.io/en/latest/>`_ makes virtual environments easier to use. It provides some very useful features:

    - it organizes the virtual environment into a single user-defined directory, so they are not scattered throughout the file system;
    - it defines commands for the easy creation, deletion, and copying of virtual environments;
    - it defines a command to activate a virtual environment using its name;
    - all commands defined by ``virtualenvwrapper`` have tab-completion for virtual environment names.

You first need to install ``virtualenvwrapper`` *outside* of a virtual environment:

.. code-block:: bash

    pip install virtualenvwrapper
    # or on ubuntu
    sudo apt install virtualenvwrapper
    # or on fedora
    sudo dnf install python-virtualenvwrapper

Then, you need to load it into your terminal session. Add the following lines in ``~/.bashrc``. They will be executed every time you open a new terminal session:

.. code-block:: bash

    # Decide where to store the virtual environments
    export WORKON_HOME=~/Envs
    # Make sure the directory exists
    mkdir -p ${WORKON_HOME}
    # Load virtualenvwrapper
    source /usr/local/bin/virtualenvwrapper.sh

Open a new terminal or run ``source ~/.bashrc`` to update your session. You can now create a virtual environment with:

.. code-block:: bash

    mkvirtualenv my-project

Regardless of your current working directory, the environment is created in ``~/Envs/`` and it is now loaded in our terminal session.

You can load your virtual environments by running ``workon my-project``, and exit them by running ``deactivate``.

Virtual environments, especially with ``virtualenvwrapper``, can do much more. For example, you can create virtual environments with different python interpreters with the ``-p`` flag. The Hitchhiker's Guide to Python has a good `tutorial <https://docs.python-guide.org/dev/virtualenvs/>`_ that gives a more in-depth explanation of virtual environments. The `virtualenvwrapper documentation <https://virtualenvwrapper.readthedocs.io/en/latest/>`_ is also a good resource to read.

On a Mac
--------

One more step is often required on macOS, because of the default number of files that a process can open simultaneously is quite low (256). To increase the number of files that can be accessed, run the following command::

    ulimit -n 4096

This sets the number of files to 4096. However, this command only applies to your currently open terminal session. To keep this high limit, add the above line to your ``~/.profile``.



.. _build-mdanalysis-develop:

-------------------
Building MDAnalysis
-------------------

Make sure that you have :ref:`cloned the repository <forking-code-repo>`  
and activated your virtual environment. First we need to install dependencies:

.. code-block:: bash

    # if using conda
    conda install cython numpy
    # if using pip
    pip install cython numpy

Ensure that you have a working C/C++ compiler (e.g. gcc or clang). You will also need Python 2.7.x or Python ≥ 3.4. We will now install MDAnalysis. 
``cd`` to the *mdanalysis/* source directory.

.. code-block:: bash

    # Build and install the MDAnalysis package
    cd package/
    python setup.py develop

    # Build and install the test suite
    cd ../testsuite/
    python setup.py develop

At this point you should be able to import MDAnalysis from your locally built version:

.. code-block:: bash

    $ python  # start an interpreter
    >>> import MDAnalysis as mda
    >>> mda.__version__
    '0.20.2-dev0'

.. _branches-in-mdanalysis:

----------------------
Branches in MDAnalysis
----------------------

There are two important branches in MDAnalysis:

    - ``master``: for production-ready code
    - ``develop``: for development code

The ``master`` branch is only for stable, production-ready code. Development code should *never* be committed to this branch. Typically, code is only committed by the release manager, when a release is ready.

The ``develop`` branch can be considered an "integration" branch for including your code into the next release. Only working, tested code should be committed to this branch. Code contributions ("features") should branch off ``develop`` rather than ``master``.


.. _create-code-branch:

Creating a branch
-----------------

The develop branch should only contain approved, tested code, so create a
feature branch for making your changes. For example, to create a branch called 
``shiny-new-feature`` from ``develop``::

    git checkout -b shiny-new-feature develop

This changes your working directory to the ``shiny-new-feature`` branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to MDAnalysis. You can have many branches with different names
and switch in between them using the ``git checkout my-branch-name`` command.

There are several special branch names that you should not use for your feature branches:

    - ``master``
    - ``develop``
    - ``release-*``
    - ``hotfix-*``


``release`` branches are used to :ref:`prepare a new production release <preparing-release>`. ``hotfix`` branches are used to :ref:`fix issues found in an already released version <preparing-hotfix>`. Both these branch types should be handled by the release manager only.

----------------
Writing new code
----------------

Code formatting in Python
-------------------------

MDAnalysis is a project with a long history and many contributors and hasn't used a consistent coding style. Since version 0.11.0 we are trying to update all the code to conform with `PEP8`_. Our strategy is to update the style every time we touch an old function and thus switch to `PEP8`_ continuously.

**Important requirements (from PEP8):**
    - keep line length to **79 characters or less**; break long lines sensibly
    - indent with **spaces** and use **4 spaces per level**
    - naming:

        - classes: `CapitalClasses` (i.e. capitalized nouns without spaces)
        - methods and functions: `underscore_methods` (lower case, with underscores for spaces)

We recommend that you use a Python Integrated Development Environment (IDE) (`PyCharm`_ and others) or external tools like `flake8`_ for code linting. For integration of external tools with emacs and vim check out `elpy`_ (emacs) and `python-mode`_ (vim).

To apply the code formatting in an automated way you can also use code formatters. External tools include `autopep8`_ and `yapf`_. Most IDEs either have their own code formatter or will work with one of the above through plugins.

Modules and dependencies
------------------------

MDAnalysis strives to keep dependencies small and lightweight. Code outside the :mod:`MDAnalysis.analysis` and :mod:`MDAnalysis.visualization` modules should only rely on the :ref:`core dependencies <core-module-dependencies>`, which are always installed. Analysis and visualization modules can use any :ref:`any package, but it is treated as optional <optional-modules>`.

Imports in the code should follow the :ref:`general-rules-for-importing`.

.. seealso::

    See :ref:`module-imports` for more information.


Developing in Cython
--------------------

The ``setup.py`` script first looks for the ``*.c`` files included in the standard MDAnalysis distribution. These are not in the GitHub repository, so ``setup.py`` will use Cython to compile extensions. ``*.pyx`` source files are used instead of ``*.c`` files. From there, ``*.pyx`` files are converted to ``*.c`` files if they are newer than the already present ``.c`` files or if the ``--force`` flag is set (i.e. ``python setup.py build --force``). End users (or devs) should not trigger the ``.pyx`` to ``.c`` conversion since ``.c`` files delivered with source packages are always up-to-date. However, devs who work on the ``.pyx`` files will automatically trigger the conversion since ``.c`` files will then be outdated. 

Place all source files for compiled shared object files into the same directory as the final shared object file.

``*.pyx`` files and cython-generated ``*.c`` files should be in the same directory as the ``*.so`` files. External dependent C/C++/Fortran libraries should be in dedicated ``src`` and ``include`` folders. See the following tree as an example::

    MDAnalysis 
        |--lib
        |   |-- _distances.so
        |   |-- distances.pyx
        |   |-- distances.c
        |-- coordinates
            |-- _dcdmodule.so
            |-- src
                |-- dcd.c
            |-- include
                |-- dcd.h


-----------------
Testing your code
-----------------

MDAnalysis takes testing seriously. All code added to MDAnalysis should have tests to ensure that it works as expected; we aim for 90% coverage. See :ref:`testing` for more on writing, running, and interpreting tests.




---------------------
Documenting your code
---------------------

Changes to the code should be reflected in the ongoing ``CHANGELOG``. Add an entry here to document your fix, enhancement, or change. In addition, add your name to the author list. If you are addressing an issue, make sure to include the issue number.



------------------------------
Adding your code to MDAnalysis
------------------------------

Committing your code
--------------------

When you are happy with a set of changes, it is time to commit. All changes in one revision should have a common theme. If you implemented two rather different things (say, one bug fix and one new feature) then split them into two commits with different messages.

Once you’ve made changes to files in your local repository, you can see them by typing::

    git status

Tell git to track files by typing::

    git add path/to/file-to-be-added.py

Doing ``git status`` again should give something like::

    # On branch shiny-new-feature
    #
    #       modified:   /relative/path/to/file-you-added.py
    #

Then commit with::

    git commit -m

This opens up a message editor. 

*Always* add a descriptive comment for your commit message (feel free to be verbose!):

    - use a short (<50 characters) subject line that summarizes the change
    - leave a blank line
    - optionally, add additional more verbose descriptions; paragraphs or bullet lists (with - or *) are good
    - manually break lines at 80 characters
    - manually indent bullet lists

.. seealso::

    See `Tim Pope's A Note About Git Commit Messages <http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html>`_ for rationale.


Pushing your code to GitHub
---------------------------

When you want your changes to appear publicly on your GitHub page, push your forked feature branch’s commits::

    git push origin shiny-new-feature

Here origin is the default name given to your remote repository on GitHub. You can see the remote repositories::

    git remote -v

If you added the upstream repository as described above you will see something like:

.. code-block:: bash

    origin	git@github.com:your-username/mdanalysis.git (fetch)
    origin	git@github.com:your-username/mdanalysis.git (push)
    upstream	git@github.com:MDAnalysis/mdanalysis.git (fetch)
    upstream	git@github.com:MDAnalysis/mdanalysis.git (push)

Now your code is on GitHub, but it is not yet a part of the MDAnalysis project. For that to happen, a pull request needs to be submitted on GitHub. 

Rebasing your code
------------------

Often the upstream MDAnalysis develop branch will be updated while you are working on your own code.
You will then need to update your own branch with the new code to avoid merge conflicts.
You need to first retrieve it and then `rebase <https://www.atlassian.com/git/tutorials/rewriting-history/git-rebase>`_
your branch so that your changes apply to the new code::

    git fetch upstream
    git rebase upstream/develop

This will replay your commits on top of the latest development code from MDAnalysis.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``git stash`` them
prior to updating.  This will effectively store your changes and they can be
reapplied after updating with ``git stash apply``. 

Once rebased, push your changes::

    git push -f origin shiny-new-feature

and `create a pull request <https://github.com/MDAnalysis/mdanalysis/pulls>`_.


Creating a pull request
-----------------------

The typical approach to adding your code to MDAnalysis is to make a `pull request`_ on GitHub:

    #. Navigate to your repository on GitHub
    #. Click on the :guilabel:`Pull Request` button
    #. You can then click on :guilabel:`Commits` and :guilabel:`Files Changed` to make sure everything looks okay one last time
    #. Write a description of your changes and follow the PR checklist

        - check that docs are updated
        - check that tests run
        - check that you've updated CHANGELOG
        - reference the issue that you address, if any

    #. Click :guilabel:`Send Pull Request`.

Your pull request is then sent to the repository maintainers. After this, the following happens:

    #. A suite of tests are run on your code with the tools `Travis`_, `Appveyor`_ and `cover?`_. If they fail, please fix your pull request by pushing updates to it.
    #. Developers will ask questions and comment in the pull request. You may be asked to make changes. 
    #. When everything looks good, a core developer will merge your code into the ``develop`` branch of MDAnalysis. Your code will be in the next release.

If you need to make changes to your code, you can do so on your local repository as you did before. Committing and pushing the changes will  update your pull request and restart the automated tests.


.. _`testing builds on Windows`: https://ci.appveyor.com/project/orbeckst/mdanalysis



Working with the documentation
==============================

MDAnalysis documentation is written in `reStructuredText <https://docutils.sourceforge.io/rst.html>`_ and built using `Sphinx`_. The
Sphinx Documentation has an excellent `introduction to reST
<http://sphinx.pocoo.org/rest.html>`__.

MDAnalysis maintains two kinds of documentation: 

    #. This user guide: a map of how MDAnalysis works, combined with tutorial-like overviews of specific topics (such as the analyses)
    
    #. The docstrings in the code itself. These are meant to provide a clear explanation of the usage of individual classes and functions. They often include technical or historical information such as in which version the function was added, or deprecation notices.

---------------------------
Working with the user guide
---------------------------

The user guide makes use of a number of Sphinx extensions to ensure that the code examples are always up-to-date. These include `nbsphinx <https://nbsphinx.readthedocs.io/en/0.5.0/>`_ and the `ipython directive <http://matplotlib.org/sampledoc/ipython_directive.html>`__.

The ``ipython`` directive lets you put code in the documentation which will be run
during the doc build. For example::

    .. ipython:: python

        x = 2
        x**3

will be rendered as::

    In [1]: x = 2

    In [2]: x**3
    Out[2]: 8

Many code examples in the docs are run during the
doc build. This approach means that code examples will always be up to date,
but it does make the doc building a bit more complex.

Forking and cloning the User Guide
----------------------------------

Go to the `MDAnalysis project page <https://github.com/MDAnalysis/UserGuide>`_ and hit the :guilabel:`Fork` button. You will
want to clone your fork to your machine:

.. code-block:: bash

    git clone https://github.com/your-user-name/UserGuide.git
    cd UserGuide
    git remote add upstream https://github.com/MDAnalysis/UserGuide

This creates the directory `UserGuide` and connects your repository to
the upstream (main project) MDAnalysis repository.


Creating a development environment
----------------------------------

:ref:`Create a new virtual environment <create-virtual-environment>` for the user guide. Install the required dependencies, and activate the ``nglview`` extension. We use ``nglview`` for visualizing molecules in Jupyter notebook tutorials.

If using conda:

.. code-block:: bash

    cd UserGuide/
    conda env create python=3.6 -f environment.yml --quiet
    conda activate mda-user-guide
    jupyter-nbextension enable nglview --py --sys-prefix


If using pip:

.. code-block:: bash

    cd UserGuide/
    pip install sphinx sphinx-sitemap nbsphinx ipython tabulate ipywidgets mdanalysis MDAnalysisTests jupyter
    pip install jupyter_contrib_nbextensions nglview
    jupyter-nbextension enable nglview --py --sys-prefix

Building the user guide
-----------------------

Navigate to the ``doc/`` directory and run ``make html``:

.. code-block:: bash

    cd doc/
    make html

The HTML output will be in ``doc/build/``, which you can open in your browser of choice. The homepage is ``doc/build/index.html``.

You can also use `sphinx-autobuild <https://pypi.org/project/sphinx-autobuild>`_ to rebuild the user guide every time you make changes to any document, including Jupyter notebooks. Install ``sphinx-autobuild``::

    pip install sphinx-autobuild

Then, run the following command in the ``doc/`` directory::

    python -m sphinx_autobuild source build

This will start a local webserver at http://localhost:8000/, which will refresh every time you save changes to a file in the user guide.

Working with Jupyter notebooks
------------------------------

One of the neat things about ``nglview`` is the ability to interact with molecules via the viewer. This ability can be preserved for the HTML pages generated from Jupyer notebooks by ``nbsphinx``, if you save the notebook widget state after execution.

------------------
Writing docstrings
------------------

The docstrings follow the **Numpy Docstring Standard**, which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.


.. _Issue Tracker: https://github.com/MDAnalysis/mdanalysis/issues
.. _`pull request`: https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests
.. _`mdnalysis-discussion`:
   http://groups.google.com/group/mdnalysis-discussion
.. _`Code of Conduct`: https://www.mdanalysis.org/pages/conduct/
.. _`Git`: http://git-scm.com/
.. _`PEP8`: https://www.python.org/dev/peps/pep-0008/
.. _`flake8`: http://flake8.readthedocs.org/en/latest/
.. _`PyCharm`: https://www.jetbrains.com/pycharm/
.. _`elpy`: https://github.com/jorgenschaefer/elpy
.. _`python-mode`: https://github.com/klen/python-mode
.. _`autopep8`: https://github.com/hhatto/autopep8
.. _`yapf`: https://github.com/google/yapf
.. _`Sphinx`: http://sphinx.pocoo.org/