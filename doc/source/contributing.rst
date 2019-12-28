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

.. note::

    Parts of this page came from the `Contributing to pandas <http://pandas.pydata.org/pandas-docs/stable/contributing.html>`_ guide.

Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

If you are brand new to MDAnalysis or open-source development, we recommend going 
through this page. If you are new to Git and version control, have a look at 
:ref:`version-control`. Get started working with the code of 
:ref:`the main MDAnalysis codebase <working-with-mdanalysis-code>` or 
:ref:`the documentation <working-with-mdanalysis-docs>`.

A brief overview of the workflow for contributing to MDAnalysis is laid out below, but it is recommended to read the whole guide for the particular section you ar interested in.

The development workflow for code or inline code documentation looks like this:

    #. :ref:`Fork the MDAnalysis repository <forking-code-repo>` from the mdanalysis account into your own account
    #. :ref:`Set up an isolated virtual environment <create-virtual-environment>` for code development
    #. :ref:`Build development versions <build-mdanalysis-develop>` of MDAnalysis and MDAnalysisTests on your computer into the virtual environment
    #. :ref:`Create a new branch off the develop branch <create-code-branch>`
    #. :ref:`Add your new feature or bug fix <writing new code>` or :ref:`add your new documentation <guidelines-for-docstrings>`
    #. :ref:`Add and run tests <testing>` (if adding to the code)
    #. :ref:`Build and view the documentation <building-code-documentation>` (if adding to the docs)
    #. :ref:`Commit and push your changes, and open a pull request. <adding-code-to-mda>`

The workflow for working with user guide looks like this:

    #. :ref:`Fork the MDAnalysis repository <forking-user-guide>` from the mdanalysis account into your own account
    #. :ref:`Set up an isolated virtual environment <create-virtual-environment-user-guide>` for your documentation
    #. :ref:`Create a new branch off the master branch <create-code-branch>`
    #. Add your new documentation.
    #. :ref:`Build and view the documentation <build-user-guide>`.
    #. :ref:`Commit and push your changes, and open a pull request <add-changes-user-guide>`.


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
<https://github.com/signup/free>`_. 

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
of MDAnalysis that do not interfere with each other, or your own stable version. 
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

Like conda, virtual environments managed with `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ allow you to use different versions of Python and Python packages for your different project. Unlike conda, virtualenv is not a general-purpose package manager; it leverages what is available on your system, and lets you install Python packages using pip.

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

Now you can install packages via pip without affecting your global environment. The packages that you install when the environment is activated will be available in terminal sessions that have the environment activated. You can deactivate the virtual environment by running::

    deactivate

The `virtualenvwrapper package <https://virtualenvwrapper.readthedocs.io/en/latest/>`_ makes virtual environments easier to use. It provides some very useful features:

    - it organises the virtual environment into a single user-defined directory, so they are not scattered throughout the file system;
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


``release`` branches are used to :ref:`prepare a new production release <preparing-release>`. ``hotfix`` branches are used to :ref:`fix issues found in an already released version <preparing-hotfix>`. Both of these branch types should be handled by the release manager only.

.. _writing-new-code:

----------------
Writing new code
----------------

Code formatting in Python
-------------------------

MDAnalysis is a project with a long history and many contributors; it hasn't used a consistent coding style. Since version 0.11.0, we are trying to update all the code to conform with `PEP8`_. Our strategy is to update the style every time we touch an old function and thus switch to `PEP8`_ continuously.

**Important requirements (from PEP8):**
    - keep line length to **79 characters or less**; break long lines sensibly
    - indent with **spaces** and use **4 spaces per level**
    - naming:

        - classes: `CapitalClasses` (i.e. capitalized nouns without spaces)
        - methods and functions: `underscore_methods` (lower case, with underscores for spaces)

We recommend that you use a Python Integrated Development Environment (IDE) (`PyCharm`_ and others) or external tools like `flake8`_ for code linting. For integration of external tools with emacs and vim, check out `elpy`_ (emacs) and `python-mode`_ (vim).

To apply the code formatting in an automated way, you can also use code formatters. External tools include `autopep8`_ and `yapf`_. Most IDEs either have their own code formatter or will work with one of the above through plugins.

Modules and dependencies
------------------------

MDAnalysis strives to keep dependencies small and lightweight. Code outside the :mod:`MDAnalysis.analysis` and :mod:`MDAnalysis.visualization` modules should only rely on the :ref:`core dependencies <core-module-dependencies>`, which are always installed. Analysis and visualization modules can use any :ref:`any package, but the package is treated as optional <optional-modules>`.

Imports in the code should follow the :ref:`general-rules-for-importing`.

.. seealso::

    See :ref:`module-imports` for more information.


Developing in Cython
--------------------

The ``setup.py`` script first looks for the ``*.c`` files included in the standard MDAnalysis distribution. These are not in the GitHub repository, so ``setup.py`` will use Cython to compile extensions. ``*.pyx`` source files are used instead of ``*.c`` files. From there, ``*.pyx`` files are converted to ``*.c`` files if they are newer than the already present ``.c`` files or if the ``--force`` flag is set (i.e. ``python setup.py build --force``). End users (or developers) should not trigger the ``.pyx`` to ``.c`` conversion, since ``.c`` files delivered with source packages are always up-to-date. However, developers who work on the ``.pyx`` files will automatically trigger the conversion since ``.c`` files will then be outdated. 

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


.. _adding-code-to-mda:

------------------------------
Adding your code to MDAnalysis
------------------------------

Committing your code
--------------------

When you are happy with a set of changes, it is time to commit. All changes in one revision should have a common theme. If you implemented two rather different things (say, one bug fix and one new feature), then split them into two commits with different messages.

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

.. _rebase-code:

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

.. _create-a-pull-request:

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

    #. A :ref:`suite of tests are run on your code <continuous-integration>` with the tools :ref:`travis`, :ref:`appveyor` and :ref:`codecov`. If they fail, please fix your pull request by pushing updates to it.
    #. Developers will ask questions and comment in the pull request. You may be asked to make changes. 
    #. When everything looks good, a core developer will merge your code into the ``develop`` branch of MDAnalysis. Your code will be in the next release.

If you need to make changes to your code, you can do so on your local repository as you did before. Committing and pushing the changes will  update your pull request and restart the automated tests.

.. _working-with-mdanalysis-docs:

Working with the documentation
==============================

MDAnalysis documentation is written in `reStructuredText <https://docutils.sourceforge.io/rst.html>`_ ("rst" or "reST") and built using `Sphinx`_. The
Sphinx Documentation has an excellent `introduction to reST
<http://sphinx.pocoo.org/rest.html>`__.

MDAnalysis maintains two kinds of documentation: 

    #. `This user guide <https://www.mdanalysis.org/UserGuide/>`__: a map of how MDAnalysis works, combined with tutorial-like overviews of specific topics (such as the analyses)
    
    #. `The documentation generated from the code itself <https://www.mdanalysis.org/docs/>`__. Largely built from code docstrings, these are meant to provide a clear explanation of the usage of individual classes and functions. They often include technical or historical information such as in which version the function was added, or deprecation notices.

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

.. _forking-user-guide:

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


.. _create-virtual-environment-user-guide:

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
    pip install -r requirements.txt
    jupyter-nbextension enable nglview --py --sys-prefix

.. _build-user-guide:

Building the user guide
-----------------------

Navigate to the ``doc/`` directory and run ``make html``:

.. code-block:: bash

    cd doc/
    make html

The HTML output will be in ``doc/build/``, which you can open in your browser of choice. The homepage is ``doc/build/index.html``.

If rebuilding the documentation becomes tedious after a while, install the :ref:`sphinx-autobuild <autobuild-sphinx>` extension. 

Saving state in Jupyter notebooks
---------------------------------

One of the neat things about ``nglview`` is the ability to interact with molecules via the viewer. This ability can be preserved for the HTML pages generated from Jupyer notebooks by ``nbsphinx``, if you save the notebook widget state after execution.

.. _add-changes-user-guide:

Adding changes to the UserGuide
-------------------------------

As with the code, :ref:`commit and push <adding-code-to-mda>` your code to GitHub. Then :ref:`create a pull request <create-a-pull-request>`. The only test run for the User Guide is: that your file compile into HTML documents without errors. As usual, a developer will review your PR and merge the code into the User Guide when it looks good.


-----------------------------------
Working with the code documentation
-----------------------------------

MDAnalysis has a lot of documentation in the Python doc strings. The docstrings follow the `Numpy Docstring Standard <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__, which is used widely
in the Scientific Python community. They are nice to read as normal text and are converted by sphinx to normal ReST through `napoleon <http://sphinxcontrib-napoleon.readthedocs.org/en/latest/index.html>`__.

This standard specifies the format of
the different sections of the docstring. See `this document
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
for a detailed explanation, or look at some of the existing functions to
extend it in a similar manner.

Note that each page of the  `online documentation <https://www.mdanalysis.org/docs/>`_ has a link to the *Source* of the page. You can look at it in order to find out how a particular page has been written in reST and copy the approach for your own documentation.

.. _building-code-documentation:

Building the documentation
--------------------------

The online documentation is generated from the pages in ``mdanalysis/package/doc/sphinx/source/documentation_pages``. The documentation for the current release are hosted at www.mdanalysis.org/docs, while the development version is at www.mdanalysis.org/mdanalysis/. 

In order to build the documentation, you must first :ref:`clone the main MDAnalysis repo <forking-code-repo>`. :ref:`Set up a virtual environment <create-virtual-environment>` in the same way as you would for the code (you can use the same environment as you do for the code). You will need to install several packages for the docs.

    .. code-block:: bash

        pip install sphinx sphinx-sitemap sphinx_rtd_theme

In addition, build the development version of MDAnalysis (if you haven't done this already)::

    python setup.py develop

Then, generate the docs with:

    .. code-block:: bash

        python setup.py build_sphinx -E

This generates and updates the files in ``doc/html``. If the above command fails with an ``ImportError``, run

    .. code-block:: bash

        python setup.py build_ext --inplace

and retry.

You will then be able to open the home page, ``doc/html/index.html``, and look through the docs. In particular, have a look at any pages that you tinkered with. It is typical to go through multiple cycles of fix, rebuild the docs, check and fix again.

If rebuilding the documentation becomes tedious after a while, install the :ref:`sphinx-autobuild <autobuild-sphinx>` extension. 

Where to write docstrings
-------------------------

When writing Python code, you should always add a docstring to each public (visible to users):

    * module
    * function
    * class
    * method
 
\When you add a new module you should include a docstring with a short sentence describing what the module does or a long document including examples and references. 

.. _guidelines-for-docstrings:

Guidelines for writing docstrings
---------------------------------

A typical function docstring looks like the following::

    def func(arg1, arg2):
        """Summary line.

        Extended description of function.

        Parameters
        ----------
        arg1 : int
            Description of `arg1`
        arg2 : str
            Description of `arg2`


        Returns
        -------
        bool
            Description of return value

        """
        return True

.. seealso::

    The `napoleon documentation <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`_ has further breakdowns of docstrings at the module, function, class, method, variable, and other levels.

* When writing reST markup, make sure that there are **at least two blank lines above** the reST after a numpy heading. Otherwise, the Sphinx/napoleon parser does not render correctly.

    .. code-block:: RST

        some more docs bla bla

        Notes
        -----
        THE NEXT TWO BLANK LINES ARE IMPORTANT.


    .. versionadded:: 0.16.0
  
* Do not use "Example" or "Examples" as a normal section heading (e.g. in module level docs): *only* use it as a `NumPy doc Section <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__. It will not be rendered as a normal section and will mess up sectioning.


* When writing multiple common names in one line, Sphinx sometimes tries to reference the first name. In that case, you have to split the names across multiple lines. See below for an example:

    .. code-block:: RST

        Parameters
        ----------
        n_atoms, n_residues : int
            numbers of atoms/residues

* We are using MathJax with sphinx so you can write LaTeX code in math tags. In blocks, the code below

    .. code-block:: rst

        #<SPACE if there is text above equation>
        .. math::
            e^{i\pi} = -1

    renders like so:

    .. math::
        e^{i\pi} = -1
    

    Math directives can also be used inline.

    .. code-block:: rst

        We make use of the identity :math:`e^{i\pi} = -1` to show...

    Note that you should *always* make doc strings with math code **raw** python strings **by prefixing them with the letter "r"**, or else you will get problems with backslashes in unexpected places. ::

        def rotate(self, R):
            r"""Apply a rotation matrix *R* to the selection's coordinates.

            :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
            :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

            .. math::

            \mathbf{x}' = \mathsf{R}\mathbf{x}
            """

    .. seealso::
    
        See `Stackoverflow: Mathjax expression in sphinx python not rendering correctly <http://stackoverflow.com/questions/16468397/mathjax-expression-in-sphinx-python-not-rendering-correclty">`_ for further discussion.


Writing docs for abstract base classes
--------------------------------------

MDAnalysis contains a number of abstract base classes, such as :class:`~MDAnalysis.analysis.base.AnalysisBase`. Developers who define new base classes, or modify existing ones, should follow these rules:

    - The *class docstring* needs to contain a list of methods that can be overwritten by inheritance from the base class. Distinguish and document methods as required or optional.
    - The class docstring should contain a minimal example for how to derive this class. This demonstrates best practices, documents ideas and intentions behind the specific choices in the API, helps to promote a unified code base, and is useful for developers as a concise summary of the API.
    - A more detailed description of methods should come in the *method docstring*, with a note specifying if the method is required or optional to overwrite.

See the documentation of :class:`MDAnalysis.analysis.base.AnalysisBase` for an example of this documentation.

Adding your documentation to MDAnalysis
---------------------------------------

As with any contribution to an MDAnalysis repository, :ref:`commit and push <adding-code-to-mda>` your documentation contributions to GitHub. If *any fixes in the restructured text* are needed, *put them in their own commit* (and do not include any generated files under `docs/html`). Try to keep all reST fixes in the one commit. ``git add FILE`` and ``git commit --amend`` is your friend when piling more and more small reST fixes onto a single "fixed reST" commit.

Then, :ref:`create a pull request <create-a-pull-request>`. All the tests in the MDAnalysis test suite will run, but only one checks that the documents compile correctly.

Viewing the documentation interactively
---------------------------------------

In the Python interpreter one can simply say::

  import MDAnalysis
  help(MDAnalysis)
  help(MDAnalysis.Universe)

In ``ipython`` one can use the question mark operator::

  MDAnalysis.Universe?

.. _autobuild-sphinx:

------------------------------------
Automatically building documentation
------------------------------------

Constantly rebuilding documentation can become tedious when you have many changes to make. Use `sphinx-autobuild <https://pypi.org/project/sphinx-autobuild>`_ to rebuild documentation every time you make changes to any document, including Jupyter notebooks. Install ``sphinx-autobuild``::

    pip install sphinx-autobuild

Then, run the following command in the ``doc/`` directory::

    python -m sphinx_autobuild source build

This will start a local webserver at http://localhost:8000/, which will refresh every time you save changes to a file in the documentation. This is helpful for both the user guide (first navigate to ``UserGuide/doc``) and the main repository documentation (navigate to ``package/doc/sphinx``).



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