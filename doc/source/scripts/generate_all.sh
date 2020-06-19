#!/usr/bin/env bash

GENPATH=$(echo $(cd $(dirname $0) && pwd -P))

python $GENPATH/gen_format_overview_classes.py
python $GENPATH/gen_selection_exporters.py
python $GENPATH/gen_standard_selections.py
python $GENPATH/gen_topology_groupmethods.py
python $GENPATH/gen_topologyattr_defaults.py
python $GENPATH/gen_topologyparser_attrs.py