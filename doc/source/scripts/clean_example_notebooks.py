#!/usr/bin/env python
"""
Use this to execute and clean up notebooks by:

    - Execute the notebook to make sure it works as intended
    - Adding Last updated
    - Numbering references
    - Adding a References section

Run the script manually with a path to the MDAnalysis repo.
"""

import os
import json
import glob
import argparse

parser = argparse.ArgumentParser(description='Clean Jupyter notebooks.')
parser.add_argument('path_to_notebooks', type=str, default='./')

def get_ipynbs(path):
    path = os.path.abspath(path)
    pattern = os.path.join(path, '**', '.ipynb')
    for file in glob.iglob(pattern, recursive=True):
        yield file

class JupyterCell:

    def __init__(self, cell_type='markdown', metadata={}, source=[]):
        self.cell_type = cell_type
        self.metadata = metadata
        self.source = source

    

class JupyterNotebook:

    def __init__(self, filename):
        with open(filename, 'r') as f:
            contents = json.load(f)

        self.cells = 

def clean_ipynb(file):
    with open(file, 'r') as f:
        contents = json.load(f)

    first_cell = 
    # insert or update Last updated
    for i, line in enumerate(contents['cells'][0]['source']):
        if 'last updated' in line.lower():
            

def reference_ipynb(contents):
    refs = []
    for cell in contents['cells']:
        pass

def read_rst()