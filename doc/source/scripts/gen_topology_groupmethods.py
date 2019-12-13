#!/usr/bin/env python
"""
Generate groupmethods.txt:

A table of transplanted methods.
"""

from collections import defaultdict
from MDAnalysis.core.groups import GroupBase
from base import TableWriter
from core import TOPOLOGY_CLS


class TransplantedMethods(TableWriter):

    headings = ['Method', 'Description', 'Requires']
    filename = 'generated/topology/groupmethods.txt'

    def _set_up_input(self):
        items = []
        for klass in TOPOLOGY_CLS:
            for name, method in klass.transplants[GroupBase]:
                items.append([name, klass, method])
        return [x[1:] for x in sorted(items)]
    
    def _method(self, klass, method):
        return self.sphinx_meth(method)
    
    def _description(self, klass, method):
        return ' '.join(method.__doc__.split('.\n')[0].split())
    
    def _requires(self, klass, method):
        return klass.attrname

if __name__ == '__main__':
    TransplantedMethods()

