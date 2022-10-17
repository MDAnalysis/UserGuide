#!/usr/bin/env bash

GENPATH=$(echo $(cd $(dirname $0) && pwd -P))

which python
python -c "import MDAnalysis as mda; print(mda.__version__)"
python -c "import MDAnalysisTests as mdat; print(mdat.__version__)"

python $GENPATH/gen_format_overview_classes.py
python $GENPATH/gen_selection_exporters.py
python $GENPATH/gen_standard_selections.py
python $GENPATH/gen_topology_groupmethods.py
python $GENPATH/gen_topologyattr_defaults.py
python $GENPATH/gen_topologyparser_attrs.py
python $GENPATH/gen_unit_tables.py
