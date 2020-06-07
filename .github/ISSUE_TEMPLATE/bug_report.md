---
name: Bug report
about: Create a report to help us fix mistakes in the user guide
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected from the user guide. Please include a link to the any sections that you reference.

**Actual behavior**
What is there, or not there, instead. Add as much detail as you can. Include (copy and paste) error messages, stack traces, and output of any code that is related.

If applicable, show us how to reproduce the failure. If you can, use trajectory files from the test data.

``` python
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD,  GRO, PDB, TPR, XTC, TRR,  PRMncdf, NCDF

u = mda.Universe(PSF, DCD)

....

```

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]
- Which version of MDAnalysis? (run `python -c "import MDAnalysis as mda; print(mda.__version__)"`)
- Which version of Python (`python -V`)?

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
