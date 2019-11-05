from __future__ import print_function
import os
import tabulate
from MDAnalysis import _TOPOLOGY_ATTRS

ignore = ('topologyattrs', 'atomattrs', 'residueattrs', 'segmentattrs', 'indices', 'resindices', 'segindices')

TOPOLOGY_CLS = sorted(set([x for x in _TOPOLOGY_ATTRS.values()
                           if not x.attrname in ignore]), 
                      key=lambda x: x.attrname)

def write_rst_table(lines, headings, filename):
    path = os.getcwd()
    if not 'scripts' in path:
        filename = os.path.join('scripts', filename)
    with open(filename, 'w') as f:
        print(tabulate.tabulate(lines, headers=headings, tablefmt='rst'), file=f)
    print('Wrote ', filename)