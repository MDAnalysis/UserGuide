#!/usr/bin/env python
"""
Use this to execute and clean up notebooks by:

    - Execute the notebook to make sure it works as intended
    - Adding Last executed line
    - Converting references
    - Adding a References section
    - Rewriting the notebook

Note: this does not preserve NGLView widgets. You will have to
re-execute those notebooks and click Widgets > Save Notebook Widget State
to keep functionality in the HTML version.

To run::

    ./clean_example_notebooks.py *.ipynb -vvv

References
==========

This script uses the sphinxcontrib-bibtex extension for references.
"""

import os
import glob
import datetime
import re
import shutil
import argparse
import logging
import sys

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import pybtex as tex
from pybtex.database import parse_file, BibliographyData
import pybtex.plugin
from pybtex.style import formatting

import MDAnalysis as mda

parser = argparse.ArgumentParser(description='Clean Jupyter notebooks.')
parser.add_argument('files', type=str, nargs='+', help='notebook files')
parser.add_argument('-v', '--verbose', action='count', default=0)


logger = logging.getLogger(__name__)


class References:
    """
    Reference manager
    """
    style_name = 'plain'
    template = '[{i}] {html}'
    ref_header = '## References'

    def __init__(self, filename='../references.bib'):
        logger.info('Reading references from {}'.format(filename))
        self.data = parse_file(filename)
        logger.info(self.data.entries.keys())

        # set up bibliography formatting
        self.style = pybtex.plugin.find_plugin('pybtex.style.formatting',
                                               self.style_name)()
        self.backend = pybtex.plugin.find_plugin('pybtex.backends',
                                                 'html')()

        # set up inline reference and doi dictionaries
        entries = self.data.entries
        self.inline = ({k: self._parse_entry(v)
                        for k, v in entries.items()})
        self.doi = ({k: self._get_doi(v)
                     for k, v in entries.items()})

        # set up key regex
        pattern = '('+')|('.join(self.data.entries.keys())+')'
        self.regex = re.compile(pattern)

    def write_bibliography(self, keys=[]):
        """
        Write ordered, numbered bibliography for each notebook.
        Output is in HTML.
        """
        bibs = [self.ref_header]
        entries = [self.data.entries[k] for k in keys]
        for i, entry in enumerate(self.style.format_entries(entries), 1):
            html = entry.text.render(self.backend)
            bibs.append(self.template.format(i=i, html=html))
        return '\n\n'.join(bibs)

    def _get_doi(self, entry):
        try:
            return 'https://doi.org/{}'.format(entry.fields['doi'])
        except KeyError:
            return '#References'

    def _parse_entry(self, entry):
        authors = entry.persons['author']
        inline = authors[0].last_names[0]
        if len(authors) > 2:
            inline += ' *et al.*'
        elif len(authors) == 2:
            inline += ' and {}'.format(authors[1].last_names[0])
        inline += ', '
        year = entry.fields['year']
        inline += year
        return inline


class JupyterCell:
    """
    Handles each Jupyter cell.
    """

    last_executed = '**Last executed:** {} with MDAnalysis {}'
    time_fmt = '%b %d, %Y'
    tag = 'a'
    close_tag = '</{}>'.format(tag)
    tagline = ('<{} data-cite="{{key}}" '
               'href="{{url}}">{{authors}}</{}>').format(tag, tag)

    @classmethod
    def as_references(cls, refs, keys=[]):
        source = refs.write_bibliography(keys)
        return cls(source=source)

    def __init__(self, cell_type='markdown', metadata={},
                 source='', **kwargs):
        self.cell_type = cell_type
        self.metadata = metadata
        self.source = source
        self.lines = source.split('\n')
        self.kwargs = kwargs

    def update_last_executed(self, now, version):
        self._update_last_executed_lines(now, version)
        self.source = '\n'.join(self.lines)

    def _update_last_executed_lines(self, now, version):
        """
        Add a **Last executed** line to the first cell.

        First tries to replace former **Last executed**.
        Then, sticks it in front of **Last updated**.
        If that is not there, then in front of **Minimum version**.
        Finally, the last option is just the last line of the first cell.
        """
        copied = self.lines[:]
        time = now.strftime(self.time_fmt)
        newline = self.last_executed.format(time, version)
        for i, line in enumerate(copied):
            if 'last executed' in line.lower():
                self.lines[i] = newline
                return
        for i, line in enumerate(copied):
            if 'last updated' in line.lower():
                self.lines.insert(i-1, '\n'+newline)
                return
        for i, line in enumerate(copied):
            if 'minimum version' in line.lower():
                self.lines.insert(i-1, '\n'+newline)
                return
        self.lines.append(newline + '\n')

    def find_reference_keys(self, refs, keys=[]):
        """
        Replace shorthand reference keys with formatted, linked
        inline references. In the Jupyter notebook these link to the
        paper DOI. In the HTML output they are displayed as Sphinx
        references.

        Track the order of the keys for a final bibliography cell.
        """
        matches = [x for x in re.split(
            refs.regex, self.source) if x is not None]
        new_source = ''

        while len(matches) > 1:
            key = matches[1]
            if key not in keys:
                keys.append(key)

            authors = refs.inline[key]
            url = refs.doi[key]

            before = matches[0]
            prev_char = before[-1]
            # shorthand
            if prev_char == '#':
                new_source += before[:-1]
            # already in an HTML tag
            elif prev_char in ('"', "'"):
                new_source += before.rsplit('<', maxsplit=1)[0]
                if len(matches) > 2:
                    matches[2] = matches[2].split(self.close_tag,
                                                  maxsplit=1)[-1]
            tag = self.tagline.format(key=key,
                                      authors=authors,
                                      url=url)
            new_source += tag
            matches.pop(0)
            matches.pop(0)

        if len(matches):
            new_source += matches.pop(0)

        self.source = new_source
        self.lines = new_source.split('\n')

    def to_dict(self):
        """
        Turn into a dictionary
        """
        cell = {
            'cell_type': self.cell_type,
            'metadata': self.metadata,
            'source': self.source
        }
        cell.update(self.kwargs)
        return cell


class JupyterNotebook:
    """
    Handles a Jupyter notebook. It immediately backs up the original to
    .original, just in case.
    """
    kernel_name = os.environ['CONDA_DEFAULT_ENV']
    version = mda.__version__

    def __init__(self, filename, refs):
        logger.info('Operating on notebook {}'.format(filename))
        with open(filename, 'r') as f:
            self.contents = nbformat.reads(f.read(), as_version=4)
        self.cells = [JupyterCell(**c) for c in self.contents['cells']]
        self.metadata = self.contents['metadata']
        self.keys = []
        self.refs = refs
        self.filename = filename
        # make backup
        split = filename.split('/')
        backup_name = '/'.join(split[:-1] + ['.'+split[-1]])
        shutil.copyfile(filename, backup_name)
        self.err = None
        self.nglview = False

    def clean(self):
        """
        Write references and execute the notebook. If an error is
        raised, save the error. If no error is raised, the notebook
        is overwritten.
        """
        self.get_references()
        try:
            self.execute_and_update()
        except Exception as err:
            self.err = err

    def get_references(self):
        """
        Get references in all the cells and build a bibliography in the
        last cell.
        """
        logger.info('   Rewriting references')
        for cell in self.cells:
            if cell.cell_type == 'markdown':
                cell.find_reference_keys(self.refs,
                                         keys=self.keys)
        if '## References' in self.cells[-1].source:
            self.cells = self.cells[:-1]

        if self.keys:
            self.cells.append(JupyterCell.as_references(self.refs,
                                                        self.keys))
        self.contents['cells'] = [c.to_dict() for c in self.cells]
        self.contents = nbformat.from_dict(self.contents)

    def execute_and_update(self):
        """
        Execute the notebook. Overwrite the notebook if
        """
        original_contents = self.contents
        logger.info('   Executing')
        # Choosing a kernel name gets me errors?
        ep = ExecutePreprocessor(timeout=600, kernel_name=self.kernel_name)
        ep.preprocess(self.contents)

        logger.info('   Updating last executed')
        # on success, replace Last executed line
        first_cell = JupyterCell(**self.contents['cells'][0])
        now = datetime.datetime.now()
        first_cell.update_last_executed(now, self.version)

        if 'nglview' in first_cell.source:
            # never mind, this doesn't successfully preserve older widget states.
            # self.contents = original_contents
            self.nglview = True

        self.contents['cells'][0] = first_cell.to_dict()
        self.contents['metadata'] = self.metadata
        self.contents = nbformat.from_dict(self.contents)

        logger.info('   Rewriting notebook')
        with open(self.filename, 'w', encoding='utf-8') as f:
            nbformat.write(self.contents, f)


def clean_all_notebooks(notebooks):
    not_backups = [n for n in notebooks if not n.split('/')[-1][0] == '.']
    if not not_backups:
        return
    refs = References()
    errs = []
    nglview = []
    for nb in not_backups:
        notebook = JupyterNotebook(nb, refs)
        notebook.clean()
        if notebook.err is not None:
            errs.append((nb, notebook.err))
        if notebook.nglview:
            nglview.append(nb)

    if nglview:
        print('Re-execute for NGLView: ')
        print('\n'.join(nglview))

    if len(errs):
        errmsgs = ['{}: {}'.format(nb, err) for nb, err in errs]
        delim = '\n'+'==='*10 + '\n'
        raise ValueError(
            'Notebooks have errors: {}'.format(delim.join(errmsgs)))


if __name__ == '__main__':
    args = parser.parse_args()
    level = max(0, 50-(args.verbose*10))
    logging.basicConfig(
        level=level,
        stream=sys.stdout,
    )
    clean_all_notebooks(args.files)
