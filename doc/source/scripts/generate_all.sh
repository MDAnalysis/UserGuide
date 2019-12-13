#!/usr/bin/env bash

python gen_format_overview_classes.py
python gen_selection_exporters.py
python gen_standard_selections.py
python gen_topology_groupmethods.py
python gen_topologyattr_defaults.py
python gen_topologyparser_attrs.py