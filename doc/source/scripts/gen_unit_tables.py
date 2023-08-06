#!/usr/bin/env python
"""
Writes unit tables
"""
import pathlib
import sys
import textwrap

import tabulate
from MDAnalysis.units import conversion_factor


def write_unit_table(filename):
    headings = ["Unit", "Conversion factor"]
    table_heading = ".. table:: {}"
    tables = []
    for data_type, items in conversion_factor.items():
        lines = sorted(list(items.items()))
        tables.append((data_type, lines))

    parent_directory = pathlib.Path(__file__).parent.parent
    parent_directory.mkdir(exist_ok=True, parents=True)
    filename = parent_directory / filename

    with filename.open("w") as f:
        f.write(".. Generated by {}\n".format(sys.argv[0]))
        for data_type, lines in tables:
            line = "\n" + "-" * len(data_type) + "\n"
            f.write(line)
            f.write(data_type.capitalize())
            f.write(line)
            f.write("\n\n")
            f.write(
                textwrap.indent(
                    tabulate.tabulate(lines, headers=headings, tablefmt="rst"), "    "
                )
            )
            f.write("\n")
    print("Wrote ", filename)


if __name__ == "__main__":
    write_unit_table("generated/units_table.txt")
