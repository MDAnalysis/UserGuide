#!/usr/bin/env bash

# -o errexit: Causes the shell to exit if any invoked command exits with a non-zero status, which signals an error. This helps catch errors and undefined variables in bash scripts.
# -o nounset: Causes the shell to exit if any variable is used in the script that is not set. This can help catch typos and other bugs.
# -o pipefail: Causes the shell to exit if any command in a pipeline fails. This prevents errors in a pipeline from being masked.
set -o errexit -o nounset -o pipefail

GENPATH=$(echo $(cd $(dirname $0) && pwd -P))

which python
python -c "import MDAnalysis as mda; print(mda.__version__)"
python -c "import MDAnalysisTests as mdat; print(mdat.__version__)"

python "${GENPATH}/gen_format_overview_classes.py"
python "${GENPATH}/gen_selection_exporters.py"
python "${GENPATH}/gen_standard_selections.py"
python "${GENPATH}/gen_topology_groupmethods.py"
python "${GENPATH}/gen_topologyattr_defaults.py"
python "${GENPATH}/gen_topologyparser_attrs.py"
python "${GENPATH}/gen_unit_tables.py"
