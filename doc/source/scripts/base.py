
from __future__ import print_function
import os
import tabulate
from collections import defaultdict

class TableWriter(object):
    """
    For writing tables with easy column switching.

    Define _set_up_input() and column methods.

    Filename relative to source.
    """

    filename = ''
    headings = []
    preprocess = []
    postprocess = []
    sort = True

    def __getattr__(self, key):
        try:
            return self.fields[key]
        except KeyError:
            return super(TableWriter, self).__getattr__(key)

    def __init__(self, *args, **kwargs):
        stem = os.getcwd().split('source')[0]
        self.path = os.path.join(stem, 'source', self.filename)
        self.fields = defaultdict(list)
        self.get_lines(*args, **kwargs)
        self.write_table()

    def _run_method(self, method, *args, **kwargs):
        sanitized = self.sanitize_name(method)
        meth = getattr(self, sanitized)
        val = meth(*args, **kwargs)
        self.fields[method].append(val)
        return val

    @staticmethod
    def sanitize_name(name):
        return '_' + name.replace(' ', '_').replace('/', '_').lower()


    def get_lines(self, *args, **kwargs):
        lines = []
        for items in self._set_up_input():
            try:
                lines.append(self.get_line(*items))
            except TypeError: # one argument
                lines.append(self.get_line(items))
        if self.sort:
            lines = sorted(lines)
        self.lines = lines

    def get_line(self, *args):
        line = []
        for p in self.preprocess:
            self._run_method(p, *args)
        for h in self.headings:
            line.append(self._run_method(h, *args))
        for p in self.postprocess:
            self._run_method(p, *args)
        return line

    def write_table(self):
        with open(self.path, 'w') as f:
            print(tabulate.tabulate(self.lines, 
                                    headers=self.headings,
                                    tablefmt='rst'),
                                    file=f)
        print('Wrote ', self.filename)


    # ==== HELPER FUNCTIONS ==== #
    
    @staticmethod
    def sphinx_class(klass, tilde=True):
        prefix = '~' if tilde else ''
        return ':class:`{}{}.{}`'.format(prefix, klass.__module__, klass.__name__)

    @staticmethod
    def sphinx_meth(meth, tilde=True):
        prefix = '~' if tilde else ''
        return ':meth:`{}{}.{}`'.format(prefix, meth.__module__, meth.__qualname__)

    @staticmethod
    def sphinx_ref(txt, label=None, add_label=True):
        if label is None:
            label = txt
        # suffix = '-label' if add_label else ''
        return ':ref:`{} <{}{}>`'.format(txt, label)
    
    @staticmethod
    def sphinx_link(txt):
        return '`{}`_'.format(txt)