#!/usr/bin/env python
"""
Writes unit tables
"""
import tabulate
import textwrap
from MDAnalysis.units import conversion_factor

def write_unit_table(filename):
    headings = ['Unit', 'Conversion factor']
    table_heading = '.. table:: {}'
    tables = []
    for data_type, items in conversion_factor.items():
        lines = sorted(list(items.items()))
        tables.append((data_type, lines))
    
    with open(filename, 'w') as f:
        for data_type, lines in tables:
            # f.write(table_heading.format(data_type.capitalize()))
            line = '\n' + '-' * len(data_type) + '\n'
            f.write(line)
            f.write(data_type.capitalize())
            f.write(line)
            f.write('\n\n')
            f.write(textwrap.indent(tabulate.tabulate(lines, 
                                    headers=headings,
                                    tablefmt='rst'),
                                    '    '))
            f.write('\n')
    print('Wrote ', filename)


if __name__ == '__main__':
    write_unit_table('generated/units_table.txt')