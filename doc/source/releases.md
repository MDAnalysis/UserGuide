# MDAnalysis Release Notes


## Release 2.7.0 of MDAnalysis

This a minor release of MDAnalysis.

This release of MDAnalysis is packaged under a [GPLv3+ license](https://www.gnu.org/licenses/gpl-3.0.en.html), additionally all contributions made from commit `44733fc214dcfdcc2b7cb3e3705258781bb491bd` onwards are made under the [LGPLv2.1+ license](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).

The minimum supported NumPy version is 1.22.3.

Supported Python versions:
  - 3.9, 3.10, 3.11, 3.12

### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.7.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.7-is-out) for more details.

#### Fixes:
* NoJump now properly handles jumps that occur on the second frame of NPT trajectories, PR #4258
* Fix charge reading from PDBQT files. PR #4283
* Fixed a case where qcprot.CalcRMSDRotationalMatrix would return a RMSD of None. PR #4273

#### Enhancements:
* Support was added for reading chainID from prmtop AMBER topologies (PR #4007)
* Added support for Python 3.12 (PR #4309, #4300, #4301, #4319, #4325, #4327, #4329)
* Added support for reading `chainID` from PDBQT files (PR #4284)
* TPR reader now sets `chainID` from `molblock` (PR #4281)
* Various improvements to the organization and performance of Major and Minor Pair analyses (PR #3735)
* C distance backend is now exposed via `libmdanalysis.pxd` (PR #4342)
* Added a GROMOS11 Reader (PR #4294)

#### Changes:
* Added `mda_xdrlib` as a core dependency to replace the now deprecated Python `xdrlib` code (PR #4271)
* ConverterBase has been moved to `MDAnalysis.converters.base` (PR #4253)
* `networkx` is now an optional dependency of MDAnalysis (PR #4331)
* `BioPython` is now an optional dependency of MDAnalysis (PR #4332)
* Results for WatsonCrickDist nucleic acids analysis are now stored in  `analysis.nucleicacids.WatsonCrickDist.results.distances` (PR #3735)

#### Deprecations:
* Importing ConverterBase from `MDAnalysis.coordinates.base` will not be possible after MDAnalysis 3.0 (PR #4253)
* Deprecation with intent of removal in MDAnalysis v3.0 of the X3DNA legacy code (PR #4333)
* Deprecation with intent of removal in MDAnalysis v3.0 of the TRZ reader and writer (PR #4335)
* Deprecation with intent of removal in MDAnalysis v3.0 of the `MDAnalysis.lib.util.which` method (PR #4340)
* The `asel` argument of the `timeseries` attribute of Readers is now deprecated in favour of the `atomgroup` argument (PR #4343)
* In `nucleicacids.WatsonCrickDist`, accepting lists of `Residue` objects was deprecated in favor of using `ResidueGroup`: using `List[Residue]` will be removed in release 3.0.0; instead use a `ResidueGroup` (PR #3735)
* In `nucleicacids.WatsonCrickDist` the result `results.pair_distances` was deprecated and will be removed in 3.0.0; use `results.distances` instead (PR #3735)

## New Contributors
* [@jennaswa](https://github.com/jennaswa) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4289
* [@Sumit112192](https://github.com/Sumit112192) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4346
* [@HeetVekariya](https://github.com/HeetVekariya) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4359
* [@JoStoe](https://github.com/JoStoe) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4292
* [@ljwoods2](https://github.com/ljwoods2) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4366

## Release 2.6.1 of MDAnalysis

This is a bugfix release of the 2.6.x version branch of MDAnalysis, it serves as an amendment to the earlier released version 2.6.0.

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.6.1/package/CHANGELOG) for more details.

### Bug fixes and changes

* Reverting the v2.6.0 behaviour, builds are now again made using the oldest supported NumPy version (NumPy 1.22.3 for Python 3.9-3.10, and 1.22.3 for Python 3.11) [PR #4261]
* Uses of numpy `in1d` have been replaced with `isin` in prepartion for NumPy 2.0 [PR #4255]
* Cython DEF statements have been replaced with compile time integer constants as DEF statements are now deprecated in Cython 3.0 [Issue #4237, PR #4246]
* Fix to element guessing code to more accurately interpret atom names split by numbers (i.e. N0A is now recognised as N rather than NA) [Issue #4167, PR #4168]
* Clarification of SurvivalProbability function documentation [Issue #4247, PR #4248]1

### New Contributors
* [@pillose](https://github.com/pillose) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4168

## Release 2.6.0 of MDAnalysis

This a minor release of MDAnalysis.

This release of MDAnalysis is packaged under a [GPLv3+ license](https://www.gnu.org/licenses/gpl-3.0.en.html), additionally all contributions made from commit 44733fc214dcfdcc2b7cb3e3705258781bb491bd onwards are made under the [LGPLv2.1+ license](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html). More details about these license changes will be provided in an upcoming blog post.

The minimum supported NumPy version has been raised to 1.22.3 as per NEP29. Please note that package builds are now made with NumPy 1.25+ which offer backwards runtime compatibility with NEP29 supported versions of NumPy.

Supported Python versions:
  - 3.9, 3.10, 3.11

### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.6.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.6-is-out) for more details.

#### Fixes:
* The -ffast-math compiler flag is no longer used by default at build time, avoiding inconsistent (although still scientifically correct) results seen in Intel MacOS systems when calling `AtomGroup.center_of_charge(..., unwrap=True). This also avoids potentially incorrect floating point results as [detailed here](https://moyix.blogspot.com/2022/09/someones-been-messing-with-my-subnormals.html). (https://github.com/MDAnalysis/mdanalysis/pull/4220)
* DATAWriter, CRD, PQR, and PDBQT files can now be correctly written to compressed files. Before this change, any attempt to write to a compressed format (gzip or bzip2) would lead to writing out an uncompressed file. (https://github.com/MDAnalysis/mdanalysis/pull/4163)
* Prevent accidental merging of bond/angle/dihedral types when they are defined as LAMMPS style string integers instead of tuples. This was leading to an incorrect number of bond/angle/dihedral types being written to lammps data files. (https://github.com/MDAnalysis/mdanalysis/pull/4003)

#### Enhancements:
* An `exclude_same` argument has been added to `InterRDF` allowing pairs of atoms that share the same residue, segment or chain to be excluded from the calculation. (https://github.com/MDAnalysis/mdanalysis/pull/4161)
* LAMMPS reader now supports the `continuous` ChainReader option. (https://github.com/MDAnalysis/mdanalysis/pull/4170)
* AtomGroup representation now returns atom indices in the same order as they are stored in the AtomGroup. (https://github.com/MDAnalysis/mdanalysis/pull/4191)

#### Changes:
* Package builds now use NumPy 1.25 or higher instead of the lowest supported NumPy version. (https://github.com/MDAnalysis/mdanalysis/pull/4198)
* As per NEP29, the minimum supported runtime version of NumPy has been increased to 1.22.3. (https://github.com/MDAnalysis/mdanalysis/pull/4160)
* The GSD package is now an optional dependency. (https://github.com/MDAnalysis/mdanalysis/pull/4174)
* The MDAnalysis package now only supports GSD versions 3.0.1 or above. (https://github.com/MDAnalysis/mdanalysis/pull/4153)
* MDAnalysis no longer officially supports 32 bit installations. (they are no longer tested in our continuous integration pipelines). Note: no code changes have been made to disable 32 bit, although it is known that new versions of most MDAnalysis core dependencies no longer release 32 bit compatible versions. (https://github.com/MDAnalysis/mdanalysis/pull/4176)
* The package license has been updated to [GPLv3+](https://www.gnu.org/licenses/gpl-3.0.en.html) to better reflect the compatibility of GPLv2+ with Apache and GPLv3 licensed codes. Additionally all new contributions from commit 44733fc214dcfdcc2b7cb3e3705258781bb491bd onwards are made under the [LGPLv2.1+ license](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html). (https://github.com/MDAnalysis/mdanalysis/pull/4219)

#### Deprecations:
* The misspelt `Boltzman_constant` entry in `MDAnalysis.units` is now deprecated in favour the correctly spelt `Boltzmann_constant`. (https://github.com/MDAnalysis/mdanalysis/pull/4230 and https://github.com/MDAnalysis/mdanalysis/pull/4214)
* `MDAnalysis.analysis.hole2` is now deprecated in favour of a new [HOLE2 MDAKit](https://www.mdanalysis.org/hole2-mdakit/). (https://github.com/MDAnalysis/mdanalysis/pull/4200)

### New Contributors
* [@MohitKumar020291](https://github.com/MohitKumar020291) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4182
* [@Shubx10](https://github.com/Shubx10) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4184
* [@ztimol](https://github.com/ztimol) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4191

## Release 2.5.0 of MDAnalysis

This a minor release of MDAnalysis.

The minimum Python version has been raised to 3.9 and NumPy to 1.21.0 as per NEP29. We also now package wheels for both linux and osx arm64 machines on PyPi.

Supported Python versions:
  - 3.9, 3.10, 3.11

### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.5.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.5-is-out) for more details.

#### Fixes:
  * Fixed an issue where transformations were not being properly applied to
    Universes with multiple trajectories (i.e. using the ChainReader)
    (Issue #3657 #4008 PR #3906)
  * Fixed an issue with the the `heavy` distance_type for `WaterBridgeAnalysis`
    where distance was not correctly assigned when more than one hydrogen
    was bonded to a heavy atom (Issue #4040, PR #4066).
  * PDB topology parser no longer fails when encountering unknown formal
    charges and instead simply does not populate attribute (Issue #4027)
  * Fixed an issue where using the `between` keyword of `HydrogenBondAnalysis` 
    led to incorrect donor-atom distances being returned (PR #4092, Issue #4091)
  * Fixed an issue where chi1_selections() ignored atom names CG1 OG OG1 SG
    and incorrectly returned `None` for amino acids CYS, ILE, SER, THR, VAL
    (Issue #4108, PR #4109)
  * Fix H5MD reader to read box vectors rather than returning `None` as the
    dimensions (Issue #4075, PR #4076)
  * Fix to allow reading NetCDF trajectories which do not have the `time`
    variable defined (Issue #4073, PR #4074)
  * Allows shape_parameter and asphericity to yield per residue quantities
    (Issue #3002, PR #3905)
  * Fix EDRReader failing when parsing single-frame EDR files (Issue #3999)
  * Add 'PairIJ Coeffs' to the list of sections in LAMMPSParser.py
    (Issue #3336)
  * PDBReader now defaults atom values for ts.data['occupancy'] to 0.0, rather
    than the previous default of 1.0. This now matches the default used when
    setting Universe Topology attributes using the first frame's information (PR #3988)

#### Enchancements:
  * ARM64 (osx and linux) wheels are now provided via PyPi (Issue #4054)
  * Addition of a new analysis class `analysis.atomicdistances.AtomicDistances`
    to provide distances between two atom groups over a trajectory.
    (Issue #3654, PR #4105)
  * Add kwarg `n_frames` to class method `empty()` in
    `MDAnalysis.core.universe`, enabling creation of a blank `Universe` with
    multiple frames (PR #4140)
  * PDBReader now populates ts.data['tempfactor'] with the tempfactor for
    each atom *for each frame*. If an entry is missing for a given atom,
    this will default to a `1.0` value. Note, this does not affect the
    topology, i.e. `AtomGroup.tempfactors` is not dynamically updated.
    (Issue #3825, PR #3988)
  * Add writing u.trajectory.ts.data['molecule_tag'] as molecule tags to
    LAMMPS data file (Issue #3548)
  * Add `progressbar_kwargs` parameter to `AnalysisBase.run` method, allowing
    to modify description, position etc of tqdm progressbars. (PR #4085)
  * Add a nojump transformation, which unwraps trajectories so that particle
    paths are continuous. (Issue #3703, PR #4031)
  * Added AtomGroup TopologyAttr to calculate gyration moments (Issue #3904,
    PR #3905)
  * Add support for TPR files produced by Gromacs 2023 (Issue #4047)
  * Add distopia distance calculation library bindings as a selectable backend
    for `calc_bonds` in `MDA.lib.distances`. (Issue #3783, PR #3914)
  * AuxReaders are now pickle-able and copy-able (Issue #1785, PR #3887)
  * Add pickling support for Atom, Residue, Segment, ResidueGroup 
    and SegmentGroup. (PR #3953)

#### Changes:
  * As per NEP29 the minimum supported Python version has been raised to 3.9
    and NumPy has been raised to 1.21 (note: in practice later versions of NumPy
    may be required depending on  your architecture, operating system, or Python
    version) (PRs #4115 and #3983).
  * Add progress bars to track the progress of `mds.EinsteinMSD` _conclude()
    methods (_conclude_simple() and _conclude_fft()) (Issue #4070, PR #4072)
  * The deprecated direct indexing and `times` from the `results` attribute of
    analysis.nucleicacids' NucPairDist and WatsonCrickDist classes has been
    removed. Please use the `results.pair_distances` and `times` attributes
    instead (Issue #3744)
  * RDKitConverter changes (part of Issue #3996):
    * moved some variables (`MONATOMIC_CATION_CHARGES `and
      `STANDARDIZATION_REACTIONS`) out of the related functions to allow users
      fine tuning them if necessary.
    * changed the sorting of heavy atoms when inferring bond orders and
      charges: previously only based on the number of unpaired electrons, now
      based on this and the number of heavy atom neighbors.
    * use RDKit's `RunReactantInPlace` for the standardization reactions, which
      should result in a significant speed improvement as we don't need to use
      bespoke code to transfer atomic properties from the non-standardized mol
      to the standardized one.

### New Contributors
* [@mglagolev](https://github.com/mglagolev) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/3959
* [@chrispfae](https://github.com/chrispfae) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4009
* [@ooprathamm](https://github.com/ooprathamm) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4010
* [@MeetB7](https://github.com/MeetB7) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4022
* [@v](https://github.com/v)-parmar made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4032
* [@MoSchaeffler](https://github.com/MoSchaeffler) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4049
* [@jandom](https://github.com/jandom) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4043
* [@xhgchen](https://github.com/xhgchen) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4037
* [@DrDomenicoMarson](https://github.com/DrDomenicoMarson) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4074
* [@AHMED](https://github.com/AHMED)-salah00 made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4059
* [@schlaicha](https://github.com/schlaicha) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4076
* [@jvermaas](https://github.com/jvermaas) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4031
* [@SophiaRuan](https://github.com/SophiaRuan) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4072
* [@marinegor](https://github.com/marinegor) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4085
* [@g2707](https://github.com/g2707) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4089
* [@DanielJamesEvans](https://github.com/DanielJamesEvans) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/4109

## Release 2.4.3 of MDAnalysis

This is a bugfix release of the 2.4.x version of MDAnalysis, it serves as an amendment to the earlier released version 2.4.2.

### Bug fixes
* Fixed DCD reading for large (>2Gb) files (Issue #4039).  This was broken for versions 2.4.0, 2.4.1 and 2.4.2.
* Fix element parsing from PSF files tests read via Parmed (Issue #4015)

## Release 2.4.2 of MDAnalysis

This is a bugfix release of the 2.4.x version of MDAnalysis, it serves as an amendment to the earlier released version 2.4.1.

### Bug fixes

* Fixed an issue where the arguments passed to np.histogramdd in `MDAnalysis.analysis.DensityAnalysis` were not compatible with the 1.24 release of NumPy (PR #3976)
* Fixed upcoming incompatibilities with NumPy 1.25 in `MDAnalysis.visualization.streamlines_3D` and `MDAnalysis.visualization.streamlines` where incorrect comparison of the truth of arrays would have led to failures (PR #3977)

## Release 2.4.1 of MDAnalysis

This is a bugfix release of the 2.4.x version of MDAnalysis, it serves as an amendment to the earlier released version 2.4.0.

### Bug fixes

* The minimum version of biopython has been raised to 1.80 for pip installs
* pytng has been added as an optional dependency 


## Release 2.4.0 of MDAnalysis

This a minor release of MDAnalysis, as per our once-every-three-months schedule.

The minimum NumPy and Python versions remain largely unchanged, however the minimum version of `biopython` has been raised to 1.80. This is also the first release to officially support Python 3.11.

Supported Python versions:
  - 3.8, 3.9, 3.10, 3.11

### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.4.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.4-is-out) for more details.

#### Fixes:

#### Enchancements:
* As part of their [outreachy project](https://www.mdanalysis.org/2022/12/15/uma_outreachy/) [@umak](https://github.com/umak) has started adding type annotations throughout the MDAnalysis codebase
* As part of their [GSoC project](https://www.mdanalysis.org/2022/11/20/bjarne_gsoc/) [@BFedder](https://github.com/BFedder) has added an auxialliary reader for EDR files (PR #3749)
* As part of their [GSoC project](https://www.mdanalysis.org/2022/12/09/Aya-gsoc-final-report/) [@aya9aladdin](https://github.com/aya9aladdin) has fixed various issues with guessing and and attribute reading. This will be followed by the introduction of a new guesser system in a future release.
* A reader for TNG files has been added by [@hmacdope](https://github.com/hmacdope), follow up on his previous [GSoC 2020](https://www.mdanalysis.org/2020/09/02/final-report-hugo/) work on creating a python library for reading TNG files (PR 3765)
* Addition of a new isolayer selection method (PR #3846)
* Various enchancements and fixes to the LAMMPS DUMP Parser (allowing box translation on reading, allowing coordinates to be unwrapped based on dump image flags, and importing of forces and velocities) (PR #3844)
* All readers now have a timeseries attribute (PR #3890)
* ReaderBase file formats now accept pathlib inputs (PR #3935)
* Added ability for hbond analysis to use types when resnames are not present (PR #3848)

#### Changes:
* The deprecated setup.py `extra_requires` `AMBER` entry has been removed in favor of `extra_formats` (PR #3810)
* Various issues with the auxilliary reader, this should not be much more robust (PR #3749)
* The Cython headers have been moved to MDAnalysis.lib.libmdanalysis (PR #3913)
* The `MDAnalysis.analysis.align.sequence_alignment` now uses Bio.Align.PairwiseAligner instead of the deprecated Bio.pairwise2 (PR #3951)

#### Deprecations:
* The MemoryReader's timeseries inclusive indexing will be changed to exclusive in version 3.0.0 (PR #3894) 
* The `sequence_aligment()` method has been deprecated and will be removed in version 3.0.0 (PR #3951)
* MDAnalysis.analysis.nucleicacids' direct indexing of selection indices to obtain pair distances results has been deprecated in favor of accessing `results.pair_distances` (PR #3958)

### New Contributors
* [@jaclark5](https://github.com/jaclark5) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/3846
* [@pgbarletta](https://github.com/pgbarletta) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/3876
* [@jfennick](https://github.com/jfennick) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/3832
* [@Hakarishirenai](https://github.com/Hakarishirenai) made their first contribution in https://github.com/MDAnalysis/mdanalysis/pull/3956



## Release 2.3.0 of MDAnalysis

This a minor release of MDAnalysis, as per our once-every-three-months schedule.

The minimum NumPy version has been raised to 1.20.0 (1.21 for macosx-arm64) in line with [NEP29](https://numpy.org/neps/nep-0029-deprecation_policy.html).


Supported python versions:
  - 3.8, 3.9, 3.10


### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.3.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.3-is-out) for more details.


#### Fixes:
 - Fixed reading error when dealing with corrupt PDB CONECT records, and an issue where MDAnalysis would write out unusable CONECT records with index>100000 (Issue #988).


#### Enhancements:
 - Formal charges are now read from PDB files and stored in a `formalcharge` attribute (PR #3755).
 - A new normalizing `norm` parameter for the `InterRDF` and `InterRDF_s` analysis methods (Issue #3687).
 - Improved Universe serialization performance (Issue #3721, PR #3710).


#### Changes:
 - To install optional packages for different file formats supported by MDAnalysis, use `pip install ./package[extra_formats]` (Issue #3701, PR #3711).


#### Deprecations:
 - The `extra_requires` target `AMBER` for `pip install ./package[AMBER]` will be removed in 2.4.0. Use `extra_formats` (Issue #3701, PR #3711).


#### CZI EOSS Performance Improvements:

A series of performance improvements to the MDAnalysis library's backend have been made as per planned work under MDAnalysis' CZI EOSS4 grant. Further details about these will be provided in a future blog post.

  - `MDAnalysis.lib.distances` now accepts `AtomGroups` as well as NumPy arrays (PR #3730). 
  - Timestep has been converted to a Cython Extension type (PR #3683).

## Release 2.2.0 of MDAnalysis

In line with NEP29, this version of MDAnalysis drops support for Python 3.7 and raises the minimum NumPy version to 1.19.0. Minimum version support has also been changed for the following packages; `networkx>=2.0`, `scipy>=1.5.0`, `gsd>=1.9.3`. Further details on MDAnalysis future support strategy and NEP29 will be released shortly.

Supported python versions:
  - 3.8, 3.9, 3.10

### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.2.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.2-is-out) for more changes and details.

#### Enhancements:
 - The `frames` argument was added to AnalysisBase-derived classes (i.e. modern analysis classes) allowing for specific frames to be defined when running an analysis. (PR #3415)
 - DL_POLY classic HISTORY files are now supported (Issue #3678)
 - Python wheels are now made available through PyPI for x86_64 architectures (Issue #1300, PR #3680)
 - Added a `center_of_charge` attribute for AtomGroups (PR #3671)
 - LinearDensity now work with UpdatingAtomGroups (Issue #2508, PR #3617)
 - Addition of a PCA transformation and an associated inverse-PCA transformation was added to the PCA analysis class (PR #3596, Issue #2703)
 - Major improvements to the RDKitConverter's accuracy (PR #3044)
   - Accuracy of 99.14% when benchmarked against ChEMBL30
   - AtomGroups containing monatomic ion charges and edge cases with nitrogen, sulfur, phosphorus and conjugated systems should now have correctly assigned bond orders and charges.
 - Addition of a new AnalysisBase derived Watson-Crick distance analysis class (PR #3611)

#### Fixes:
 - Fixed issues where calling the `copy` method of Readers did not preserve optional arguments (Issue #3664, PR #3685)
 - Fixed several issues where iterating trajectories had undefined behaviour
   - Iterating (not in memory) SingleFrame readers now reset modified trajectory attributes (Issue #3423)
   - Iterating using defined indices did not rewind the trajectory (Issue #3416)
 - Fixed issues with competing processes writing to an XTC offset file leading to offset corruption (Issue #1988, PR #3375)
 - Fixed issue preventing OpenMMTopologyParsers from parsing systems with missing elements (Issue #3317, PR #3511)
 - Fixed issue with encore.covariance.covariance_matrix not working when providing an external reference (Issue #3539, PR #3621)
 - Fixed issue with broken code paths for "residues" and "segment" groupings for LinearDensity (Issue #3571, PR #3572)
 - Improved the flexibility of MOL2 reading, allowing for optional columns (`subst_id`, `subst_name` and `charge`) not to be provided (Issue #3385, PR #3598)
 - Fixed several issues related to converting AtomGroups to RDKit molecules (PR #3044):
   - Atoms are now in the same order
   - `atom.GetMonomerInfor().GetName()` now follows the guidelines for PDB files
   - Using `NoImplicit=False` no longer throws a `SanitizationError`
 - Fixed issues with incorrect reading of triclinic boxes from DUMP files (Issue #3386, PR #3403)
 - Fixed issue with the BAT method modifying input coordinate data (Issue #3501)

#### Changes:
 - The number of matches allowed when doing a smarts selection has been increased from the default
   1000 to max(1000, n_atoms * 10), an additional set of `smarts_kwargs` can now also be passed
   to override this behaviour (Issue #3469, PR #3470)
 - The `fasteners` package is now a core dependency (PR #3375)
 - LinearDensity now saves the histogram bin edges for easier plotting as `hist_bin_edges for
   each dimension in the results dictionary (Issue #2508, PR #3617)
 - ContactAnalysis now accepts AtomGroups (Issue #2666, PR #3565)

#### Deprecations:
 - The following results attribute for LinearDensity are now deprecated: (Issue #2508, PR #3617)
   - `pos` is now `mass_density`
   - `char` is now `charge_density`
   - `std` entries are now `stddev`

### Known test failures:
 - Windows builds
    * In some cases users may get permission errors with tests involving symlinks. This should not impact the behaviour of MDAnalysis but may impact the creation of temporary files when using HOLE2 (see: https://github.com/MDAnalysis/mdanalysis/issues/3556).

## Release 2.1.0 of MDAnalysis

In line with ongoing attempts to align with NEP29, this version of MDAnalysis drops support for Python 3.6 and raises the minimum NumPy version to 1.18.0.

Please note that at time of release whilst all the MDAnalysis core functionality supports Python 3.10, some optional modules do not due to a lack of support by dependencies which they require. We hope that this support will gradually be added as more of these dependencies release new versions compatible with Python 3.10.

Supported python versions:
  - 3.7, 3.8, 3.9, 3.10

### Major changes:

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.1.0/package/CHANGELOG) and our [release blog post](https://www.mdanalysis.org/blog/#mdanalysis-2.1-is-out) for more changes and details.

#### Enhancements:
 - Addition of a new dielectric analysis module (PR #2118)
 - The TPR parser now supports reading files from GROMACS 2022 (PR #3514)
 - The H5MDReader can now load trajectories without a topology (PR #3466)
 - Custom compiler flags can be used when building MDAnalysis from source (PR #3429)
 - The RDKit reader now supports parsing R/S chirality (PR #3445)
 - A new method to apply the minimum image convention to a collection of vectors, `minimize_vectors`, has been introduced (PR #3472)

#### Fixes:
 - Fixed various integer overflow issues in the distance calculation backend of MDAnalysis which would prevent calculations on large systems (Issues #3183, #3512).
 - Fixed issues with the creation of VMD surfaces in HOLE2 when using a non-contiguous start/stop/step.
 - Fixes reading of charges with the ITPParser (Issue #3419).
 - Fixed issue with the creation of a Universe from a custom object which only provides a topology (Issue #3443).
 - Fixed issue with accessing newly created values added via `add_Segment` or `add_Residue` (Issue #3437).

#### Changes:
 - `packaging` is now a core dependency of MDAnalysis.
 - Indexing a Group (AtomGroup, ResidueGroup, SegmentGroup) with `None` now raises a `TypeError`. Prior to this indexing by `None` would incorrectly return the whole Group but claim to have a length of 1 atom (Issue #3092).
 - The TRZReader now defaults to a `dt` value of 1.0 ps instead of the previous 0.0 ps (Issue #3257).

#### Deprecations:
 - The `pbc` keyword argument for various Group methods has been deprecated in favor of `wrap`. The deprecated keyword will be removed in version 3.0.0 (Issue #1760).

### Known test failures:
 - `pytest-xdist` and more than 4 workers
   * Under these conditions a test related to logging for HydrogenBondAnalysis can fail. This is not thought to impact the validity of MDAnalysis. See here for more details: https://github.com/MDAnalysis/mdanalysis/issues/3543
 - Windows builds
    * In some cases users may get permission errors with tests involving symlinks. This should not impact the behaviour of MDAnalysis but may impact the creation of temporary files when using HOLE2 (see: https://github.com/MDAnalysis/mdanalysis/issues/3556).

## Release 2.0.0 of MDAnalysis

This is the first version of MDAnalysis to solely support python 3.6+

Supported python versions:
  - 3.6, 3.7, 3.8, 3.9

Please note that starting with the next minor version, MDAnalysis will be following [NEP29](https://numpy.org/neps/nep-0029-deprecation_policy.html).

### Notes:
* This is a major release and introduces major advertised API breaks. Caution is advised when upgrading to 2.0.0.

### Major changes:

#### Enhancements:
- LAMMPSDumpReader can now read coordinates in all different LAMMPS coordinate conventions (Issue #3358)
- New `Results` class for storing analysis results (Issue #3115)
- New OpenMM coordinate and topology converters (Issue #2863, PR #2917)
- New `intra_bonds`,`intra_angles`, `intra_dihedrals`, etc... methods to return connections involve atoms within AtomGroups instead of including atoms outside of it (Issue #1264, #2821, PR #3200)
- Support for Groamcs 2021 TPR files (Issue #3180)
- Adds preliminary support for ppc64le and aarch64 [ARM] (Issue #3127, PR #2956 #3149)
- New selection operators (Issue #3054, PR #2927)
- New refactor of helix analysis class as `analysis.helix_analysis` (Issue #2452)
- New converter between RDKit molecules and MDAnalysis AtomGroup objects (Issue #2468). Also includes `from_smiles` Universe generator method, and the aromatic and smarts selection.
- New analysis method for calculating Mean Squared Dsiplacements (Issue #2438)
- New converter between Cartesian and Bond-Angle-Torsion coordinates (PR #2668)
- Universes and readers can now be `pickled` paving the way to easier parallel analyses (Issue #2723)
- New H5MDReader and H5MDWriter (Issue #762, #2866)

#### Fixes:
- Fixes an issue where `select_atom`, AtomGroup.unique, ResidueGroup.unique, and SegmentGroup.unique did not sort the output atoms (Issues #3364 #2977)
- GRO files now only support unit cells defined with 3 or 9 entries (Issue #3305)
- Fixes the sometimes wrong sorting of atoms into fragments when unwrapping (Issue #3352)
- Fixes issue when atttempting to use/pass mean positions to PCA analysis (Issue #2728)
- Fixes support for DL_POLY HISTORY files that contain cell information even if there are no periodic boundary conditions (Issue #3314)
- Fixes issue with WaterBridgeAnalysis double counting waters (Issue #3119)
- PDBWriter will use chainID instead of segID (Issue #3144)
- PDBParser and PDBWriter now assign and use the element attribute (Issues #3030 #2422)
- AtomGroup.center now works correctly for compounds + unwrapping (Issue #2984)
- Documents and fixes the `density` keyword for `rdf.InterRDF_s` (Isuse #2811)
- Fixed Janin analysis residue filtering, including CYSH (Issue #2898)

#### Changes:
- New converter API for all MDAnalysis converters under MDAnalysis.converters
- Timestep now stores information in 'C' memory layout instead of the previous 'F' default (PR #1738)
- `hbonds.hbond_analysis` has been remove din favour of `hydrogenbonds.hbond_analysis` (Issues #2739, #2746)
- TPRParser now loads TPR files with `tpr_resid_from_one=True` by deafult, which starts TPR resid indexing from 1 (instead of 0 as in previous MDAnalysis versions) (Issue #2364, PR #3152)
- `analysis.hole` has now been removed in favour of `analysis.hole2.hole` (Issue #2739)
- `Writer.write(Timestep)` and `Writer.write_next_timestep` have been removed. Please use `write()` instead (Issue #2739)
- Removes deprecated density_from_Universe, density_from_PDB, Bfactor2RMSF, and notwithin_coordinates_factory from MDAnalysis.analysis.density (Issue #2739)
- Changes the minimum numpy supported version to 1.16.0 (Issue #2827)
- Removes deprecated `waterdynamics.HydrogenBondLifetimes` (PR #2842)
- `hbonds.WaterBridgeAnalysis` has been moved to `hydrogenbonds.WaterBridgeAnalysis` (Issue #2739 PR #2913)

#### Deprecations:
- The `bfactors` attribute is now aliased to `tempfactors` and will be removed in 3.0.0 (Issue #1901)
- `WaterBridgeAnalysis.generate_table()` now returns table information, with the `table` attribute being deprecated
- Various analysis result attributes which are now stored in `Results` will be deprecated in 3.0.0 (Issue #3261)
- In 3.0.0 the ParmEd classes will only be accessible from the `MDAnalysis.converters` module
- In 2.1.0 the TRZReader will default to a dt of 1.0 ps when failing to read it from the input TRZ trajectory

See the [CHANGELOG](https://github.com/MDAnalysis/mdanalysis/blob/release-2.0.0/package/CHANGELOG) for more changes and details.

### Known issues:
 - Windows builds
    * For some compilers (seen on MVC v.19xx), differences in floating point precision leads to PBC wrapping differing from expected outcomes. This leads to failures in the `MDAnalysisTests.lib.test_augment` tests. To our knowledge this does not significantly affect results (as all other tests pass). We will aim to fix this in version 2.1.0.

