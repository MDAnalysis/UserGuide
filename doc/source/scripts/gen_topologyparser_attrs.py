#!/usr/bin/env python
"""
Generate three tables:
    - connectivityattrs.txt: A table of supported formats for bonds, angles, dihedrals, impropers
    - topologyattrs.txt: A table of supported formats for non-connectivity attributes.
    - topology_parsers.txt: A table of all formats and the attributes they read and guess

This script imports the testsuite, which tests these.
"""
import os
import importlib
import pkgutil
import inspect
import sys
from collections import defaultdict
from core import DESCRIPTIONS, NON_CORE_ATTRS
from base import TableWriter
from _pytest.outcomes import Skipped

import MDAnalysisTests.topology
from MDAnalysisTests.topology.base import mandatory_attrs

MANDATORY_ATTRS = set(mandatory_attrs) 

def collect_classes():
    """Collect test classes of each format. Gross hacking"""
    
    parser_attrs = {}

    # iterate over all files in MDAnalysisTests.topology
    for _, name, __ in pkgutil.iter_modules(MDAnalysisTests.topology.__path__):
        try:
            module = importlib.import_module(f"MDAnalysisTests.topology.{name}")
        except Skipped:
            continue
        for name, obj in inspect.getmembers(module):
            if name.startswith("Test") and inspect.isclass(obj):
                try:
                    parser = getattr(obj, "parser")
                except AttributeError:
                    continue
                else:
                    if parser not in parser_attrs:
                        parser_attrs[parser] = {"expected_attrs": set(),
                                                "guessed_attrs": set(),}
                
                for attr in parser_attrs[parser]:
                    try:
                        parser_attrs[parser][attr] |= set(getattr(obj, attr))
                    except AttributeError:
                        pass

    for val in parser_attrs.values():
        val["expected_attrs"] = val["expected_attrs"] - MANDATORY_ATTRS
    
    return parser_attrs

class TopologyParsers(TableWriter):
    headings = ['Format', 'Description', 'Attributes read', 'Attributes guessed']
    preprocess = ['keys']
    filename = 'formats/topology_parsers.txt'
    include_table = 'Table of supported topology parsers and the attributes read'
    sort = True

    def __init__(self):
        self.attrs = defaultdict(set)
        super(TopologyParsers, self).__init__()

    def _set_up_input(self):
        parser_attrs = collect_classes()
        inp = [(x, y["expected_attrs"], y["guessed_attrs"]) for x, y in parser_attrs.items()]
        return inp

    def get_line(self, parser, expected, guessed):
        line = super(TopologyParsers, self).get_line(parser, expected, guessed)
        for a in expected|guessed:
            self.attrs[a].add(self.fields['Format'][-1])
        return line
    
    def _keys(self, parser, *args):
        f = parser.format
        if isinstance(f, (list, tuple)):
            key = f[0]
            label = ', '.join(f)
        else:
            key = label = f
        return (key, label)

    
    def _description(self, *args):
        key, label = self.keys[-1]
        return DESCRIPTIONS[key]

    
    def _format(self, *args):
        key, label = self.keys[-1]
        return self.sphinx_ref(label, key, suffix='-format')
    
    def _attributes_read(self, parser, expected, guessed):
        vals = sorted(expected - guessed)
        return ', '.join(vals)
    
    def _attributes_guessed(self, parser, expected, guessed):
        return ', '.join(sorted(guessed))
    

class TopologyAttrs(TableWriter):

    headings = ('Atom', 'AtomGroup', 'Description', 'Supported formats')
    filename = 'generated/topology/topologyattrs.txt'

    def __init__(self, attrs):
        self.attrs = attrs
        super(TopologyAttrs, self).__init__()

    def _set_up_input(self):
        return sorted([x, *y] for x, y in NON_CORE_ATTRS.items() if x not in MANDATORY_ATTRS)
    
    def _atom(self, name, singular, *args):
        return singular
    
    def _atomgroup(self, name, *args):
        return name
    
    def _description(self, name, singular, description):
        return description
    
    def _supported_formats(self, name, singular, description):
        return ', '.join(sorted(self.attrs[name]))

class ConnectivityAttrs(TopologyAttrs):
    headings = ('Atom', 'AtomGroup', 'Supported formats')
    filename = 'generated/topology/connectivityattrs.txt'

    def _set_up_input(self):
        inp = [[x]*3 for x in 'bonds angles dihedrals impropers'.split()]
        return inp


if __name__ == '__main__':
    top = TopologyParsers()
    TopologyAttrs(top.attrs)
    ConnectivityAttrs(top.attrs)