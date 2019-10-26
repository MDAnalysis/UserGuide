#!/usr/bin/env python
from __future__ import print_function
import os
import tabulate
from MDAnalysis.core import selection as sel

SELECTIONS = {
    'protein_selections': (sel.ProteinSelection, 'prot_res', True),
    'protein_backbone_selections': (sel.BackboneSelection, 'bb_atoms'),
    'nucleic_selections': (sel.NucleicSelection, 'nucl_res'),
    'nucleic_backbone_selections': (sel.NucleicBackboneSelection, 'bb_atoms'),
    'base_selections': (sel.BaseSelection, 'base_atoms'),
    'sugar_selections': (sel.NucleicSugarSelection, 'sug_atoms')
}

def chunk_list(lst, n=8):
    return [lst[i:i+n] for i in range(0, len(lst), n)]

def write_tabulated_selection(name, cls, attr, sort=False, n=8):
    selected = getattr(cls, attr)
    if sort:
        selected = sorted(selected)
    
    table = chunk_list(list(selected), n=n)
    filename = os.path.join('scripts', name+'.txt')
    with open(filename, 'w') as f:
        print(tabulate.tabulate(table, tablefmt='grid'), file=f)


if __name__ == '__main__':
    for name, args in SELECTIONS.items():
        write_tabulated_selection(name, *args)