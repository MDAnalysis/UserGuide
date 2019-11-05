#!/usr/bin/env python
"""
Generate topology_groupmethods.txt:

A table of transplanted methods.
"""

from collections import defaultdict
from MDAnalysis.core.groups import GroupBase
from base import TOPOLOGY_CLS, write_rst_table

HEADINGS = ['Method', 'Description', 'Requires']

def get_lines():
    lines = []
    for klass in TOPOLOGY_CLS:
        gb = klass.transplants[GroupBase]
        for name, method in gb:
            desc = ' '.join(method.__doc__.split('.\n')[0].split())
            mstr = '.'.join([method.__module__, klass.__name__, method.__name__])
            mstr = ':meth:`~' + mstr + '`'
            lines.append((name, mstr, desc, klass.attrname))
    return [x[1:] for x in sorted(lines)]

if __name__ == '__main__':
    lines = get_lines()
    write_rst_table(lines, HEADINGS, 'topology_groupmethods.txt')


