.. -*- coding: utf-8 -*-
.. _preparing-release:

===================
Preparing a release
===================

Rules for release branches:

    - Branch from ``develop``
    - Naming convention: ``package-*`` where ``*`` is a version number

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

Summary of tasks
----------------

  * Declare a `feature freeze` on `develop` via discord and/or the developer mailing list
  * Finalize the ``CHANGELOG`` file with the date of release
  * Increment the version across ``MDAnalysis`` and ``MDAnalysisTests`` (4 places)
  * `Merge` changes into `develop`
  * Create a tag from `develop` named `release-<version_number>`
  * Create a `package-<version_number>` branch from `develop`
  * Add rebuilt `C / C++` files to `package-<version_number>`
  * Create a tag from `package-<version_number>`
  * Check automated testing of source distribution and wheel generation and upload
  * Create a new release from newly created tag
  * Check that deployment actions have adequately pushed dist and wheels to PyPi
  * Manually upload Cirrus CI wheels (temporary)
  * Update `conda-forge` packages
  * Create a blog post outlining the release
  * Increment develop branch files ready for the next version
  * Clean up old dev docs builds


Getting the develop branch ready for a release
----------------------------------------------

#. Declare feature freeze on ``develop`` via `discord` and the `developer mailing list`_

#. Create a pre-release feature branch from ``develop``

#. Finalise the ``CHANGELOG`` with the release number and date. Ensure that the ``CHANGELOG`` summarizes important changes and includes all authors that contributed to this release.

#. Make sure the version number matches the release version. The following files need to be updated: ``package/MDAnalysis/version.py``, ``package/setup.py``, ``testsuite/MDAnalysisTests/__init__.py``, and ``testsuite/setup.py``.

#. Create a pull request against ``develop`` from this branch.


Packaging the release
---------------------

#. Create a new tag from ``develop`` named ``release-<version_number>`` where ``<version_number>`` is the release version number (this tag contains a snapshot of the Python source files as they were when the release was created):

    .. code-block:: bash

        git tag -m "release 0.7.5 of MDAnalysis and MDAnalysisTests" release-0.7.5
        git push --tags origin

#. Create a ``package-<version_number>`` branch from ``develop``. This branch is automatically protected, so you will also need to create a separate branch to create commits via PR against ``package-<version_number>``.

#. Generate new C/C++ files and commit them to the ``package-<version_number>`` branch via PR. We recommend generate the C/C++ files by building a *source distribution* tarball so you can check at the same time that there are no issues with the distribution creation process:

    .. code-block:: bash

        # MDAnalysis
        cd package/
        pipx run build --sdist

#. Once committed, create a new tag based on ``package-<version_number>`` (this tag will contain a record of all the files as they were deployed to users for that version):

    .. code-block:: bash

        git tag -m "package 0.7.5 of MDAnalysis and MDAnalysisTests" package-0.7.5
        git push --tags origin

#. Upon creation of the new ``package-*`` tag, the `deploy github action`_ workflow will be automatically triggered to create source/wheels, upload them to testpypi, re-download them and run tests.

#. If all the tests come back green, you are good to go for a full release.

    #. If tests fail you will need to work out the cause of the failure.

        #. A temporary github actions failure

            Re-run the action and wait for the tests to complete

        #. An issue with the source code.

            #. Delete the current ``package-*`` branch, and the newly created tags

            #. Add the new changes to ``develop`` and restart the release process.

            #. If the code had successfully uploaded to testpypi and failed later, you will need to create a test ``package-*`` tag which contains a different release number of in the source code (bumpy by a minor release or add a ``-beta`` modifier). Note: if the code had not successfully uploaded you can just continue the release process as normal.

            #. If CI comes back green then delete the test tag, and create a normal ``package-*`` tag with the correct version number.

            #. The github action will fail, but this is ok since we tested it with the test tag above.


Completing the release
----------------------

If everything works, you can now complete the release by:

#. Creating a release on GitHub based on the newly created ``package-<version_number>`` tag.

#. Make sure you include relevant release notes, including any known issues and highlights for the release.

#. Once published, the `deploy github action`_ will be triggered which will upload the source distributions and wheels to PyPI.

    #. If the `deploy github action`_ fails and no files have been uploaded, then restart the action.

    #. If the action fails and some files have been uploaded, then you will not be able to re-upload to PyPI. At this point you will need to yank the release from PyPI and create a new minor version and re-deploy it.


Manually upload Cirrus CI wheels (temporary)
--------------------------------------------


.. todo:: actually add some examples & links of how to do this.`


Unfortunately the deployment of Cirrus CI generated wheels (for `osx-arm64` and `linux-aarch64`) does not get properly triggered by a release.

#. Go to the ``package-*`` tag triggered Cirrus CI run results and download the generated wheels.

#. Upload them to PyPi using ``twine``.


Update `conda-forge` packages
-----------------------------

On push to PyPI, the conda-forge bot should automatically pick up the presense of a new version and create a pull request on the `MDAnalysis feedstock`_ and the `MDAnalysisTests feedstock`_. You will need to merge the MDAnalysis feedstock followed by the MDAnalysisTests feedstock in order for the new package to appear on conda-forge.

To do this you will need to:

#. Update the ``meta.yaml`` files as necessary, especially bumping up the python and dependency minimum versions as necessary.

#. If NumPy pins differ from those conda-forge uses, you will need to update the ``conda_build_config.yaml`` accordingly.

#. Ask the conda-forge bot to re-render, check that CI returns green, approve and merge the pull request.


Create a blog post outlining the release
----------------------------------------

Create a blog post outlining the release notes and publicize it on the mailing list / discord / twitter/ etc...!


Increment develop branch files ready for the next version
---------------------------------------------------------

Once the release is completed you can go ahead and update the ``develop`` branch so that it is ready for the next round of development.

#. Update the 4 version file locations with the ``-dev0`` appended version of the next release.

#. Update the CHANGELOG with a new entry for the next release.

#. Once these changes are merged into the ``develop`` branch, message the developers on discord and the mailing list letting them know that the feature freeze is over.


Clean up old developer builds of the documentation
--------------------------------------------------

Whilst new docs are automatically deployed on a release, old developer builds (appended with ``-dev``) are not automatically cleaned up. To avoid causing large amounts of files being uploaded to GitHub Pages, we need to delete these old developer builds manually. To do this switch to the ``gh-pages`` branch, delete these old files, and push the change directly.


.. _`developer mailing list`: https://groups.google.com/forum/#!forum/mdnalysis-devel
.. _`deploy github action`: https://github.com/MDAnalysis/mdanalysis/tree/develop/.github/workflows/deploy.yaml
.. _`MDAnalysis feedstock`: https://github.com/conda-forge/mdanalysis-feedstock
.. _`MDAnalysisTests feedstock`: https://github.com/conda-forge/mdanalysistests-feedstock
.. _`stable branch of the docs page`: https://docs.mdanalysis.org/stable
