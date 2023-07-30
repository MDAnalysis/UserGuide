.. -*- coding: utf-8 -*-
.. _preparing-release:

===================
Preparing a release
===================

Rules for a release branch:

    - May branch from: ``develop``
    - Must be merged into: ``master`` (and ``develop`` if needed)
    - Naming convention: ``release-*`` where ``*`` is a version number

Release policy and release numbering
====================================

We use a **MAJOR.MINOR.PATCH** scheme to label releases. We adhere to the idea of `semantic versioning <http://semver.org/>`_ (semantic versioning was introduced with release 0.9, see `Issue 200`_): Given a version number **MAJOR.MINOR.PATCH**, we increment the:

  * **MAJOR** version when we make **incompatible API changes**,
  * **MINOR** version when we **add functionality** in a backwards-compatible manner, and
  * **PATCH** version when we make backwards-compatible **bug fixes**.

However, as long as the **MAJOR** number is **0** (i.e. the API has not stabilized), even **MINOR** increases *may* introduce incompatible API changes. As soon as we have a 1.0.0 release, the public API can only be changed in a backward-incompatible manner with an increase in MAJOR version.

Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

The `CHANGELOG <https://github.com/MDAnalysis/mdanalysis/blob/develop/package/CHANGELOG>`_ lists important changes for each release.

*MAJOR*, *MINOR* and *PATCH* numbers are integers that increase monotonically.

The **release number** is set in `setup.py <https://github.com/MDAnalysis/mdanalysis/blob/develop/package/setup.py>`_ *and* in ``MDAnalysis.__version__`` (MDAnalysis/version.py), e.g. ::

    RELEASE = '0.7.5'


**While the code is in development** (i.e. whenever we are not preparing a release!) the release number gets the suffix **-dev0**, e.g. ::

    RELEASE = '0.7.6-dev0'

so that people using the :ref:`develop branch <branches-in-mdanalysis>` from the source repository can immediately see that it is not a final release. For example, "0.7.6-dev0" is the state *before* the 0.7.6 release.

.. _`Issue 200`: https://github.com/MDAnalysis/mdanalysis/issues/200

Typical workflow for preparing a release
========================================

#. Declare feature freeze on ``develop`` via the `developer mailing list`_

#. Create a release branch from ``develop``::

    git checkout -b release-0.7.6 develop

#. Finalise the ``CHANGELOG`` with the release number and date. Summarize important changes and add all authors that contributed to this release.

#. Make sure the version number is right. The following files need to be updated: ``maintainer/conda/MDAnalysis/meta.yaml``, ``package/MDAnalysis/version.py``, ``package/setup.py``, ``testsuite/MDAnalysisTests/__init__.py``, and ``testsuite/setup.py``.

#. Check that the documentation is up-to-date and tests pass.

#. Commit the finalized state::

    git commit -m "release 0.7.6 ready"

#. Create a pull request against ``develop`` from this branch.

#. Create a new tag from ``develop`` named ``release-<version_number>`` where ``<version_number>`` is the release version number (this tag contains a snapshot of the Python source files as they were when the release was created):

    .. code-block: bash

        git tag -m "release 0.7.5 of MDAnalysis and MDAnalysisTests" release-0.7.5
        git push --tags origin

#. Create a ``package-<version_number>`` branch from ``develop``. This branch is automatically protected, so you will also need to create a separate branch to create commits via PR against ``package-<version_number>``.

#. Generate new C/C++ files and commit them to the ``package-<version_number>`` branch via PR. We recommend generate the C/C++ files by building a *source distribution* tarball so you can check at the same time that there are no issues with the distribution creation process:

    .. code-block:: bash

        # MDAnalysis
        cd package/
        python setup.py sdist

#. Once committed, create a new tag based on ``package-<version_number>`` (this tag will contain a record of all the files as they were deployed to users for that version):

    .. code-block:: bash

        git tag -m "package 0.7.5 of MDAnalysis and MDAnalysisTests" package-0.7.5
        git push --tags origin

#. Upon creation of the new ``package-*`` tag, the `developer mailing list`_, will be triggered to create source/wheels, upload them to testpypi, and re-download them and run tests.

#. If all the tests come back green, you are good to go for a full release.

    #. If tests fail you will need to work out the cause of the failure.

        #. A temporary github actions failure

            Re-run the action and wait for the tests to complete

        #. An issue with the source code.

            #. Delete the current ``package-*`` branch, and the newly created tags

            #. Add the new changes to ``develop`` and restart the release process.

            #. Upon getting to the ``package-*`` tag creation, create a test tag with a different version number (bump by a minor release or add a ``-beta`` modifier), so as to trigger the release and tests against testpypi.

            #. Delete the test tag, create a normal ``package-<version_number>`` tag with the correct version number.

            #. The githuba action will fail, but this is ok since we tested it with the test tag above.

#. If everything works, you can now complete the release by createing a release on Github based on the newly created ``package-<version_number`` tag.

#. This will re-trigger the `deploy github action`_ and upload the source distributions / wheels to PyPI.

    #. If the action fails as this point and no files have been deployed, then restart the action.

    #. If the action fails and some files have been deployed, then you will not be able to re-upload to PyPI. At this point you will need to yank the release from PyPI and create a new minor version and re-deploy it.

#. Update the release on conda-forge

    #. On push to PyPI, the conda-forge bot should automatically pick up the presense of a new version and create a pull request on the `MDAnalysis feedstock`_

    #. Update the ``meta.yaml`` file as necessary, especially bumping up the python and dependency minimum versions as necessary.

    #. Ask the conda-forge admin to re-render, check that CI returns green, approve and merge the pull request.

#. Documentation will be built automatically and versioned. Check that these have been created appropriately on the `stable branch of the docs page`_.

#. Create a blog post outlining the release notes and publicize it on the mailing list / discord / twitter/ etc...!

.. _`developer mailing list`: https://groups.google.com/forum/#!forum/mdnalysis-devel
.. _`deploy github action`: https://github.com/MDAnalysis/mdanalysis/tree/develop/.github/workflows/deploy.yaml
.. _`MDAnalysis feedstock`: https://github.com/conda-forge/mdanalysis-feedstock
.. _`stable branch of the docs page`: https://docs.mdanalysis.org/stable
