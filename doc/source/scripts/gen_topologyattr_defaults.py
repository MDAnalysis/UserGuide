"""
Generate topology_defaults.txt:

A table of whether TopologyAttrs are atomwise, residuewise, or segmentwise, and their defaults
"""

from base import TableWriter
from core import TOPOLOGY_CLS
from MDAnalysis.core.topologyattrs import (
    AtomAttr,
    ResidueAttr,
    SegmentAttr,
    TopologyAttr,
)

DEFAULTS = {
    "resids": "continuous sequence from 1 to n_residues",
    "resnums": "continuous sequence from 1 to n_residues",
    "ids": "continuous sequence from 1 to n_atoms",
}


class TopologyDefaults:
    def __init__(self) -> None:
        def _atom(klass: TopologyAttr) -> str:
            return klass.attrname  # type: ignore

        def _atomgroup(klass: TopologyAttr) -> str:
            return klass.singular  # type: ignore

        def _default(klass: TopologyAttr) -> str:
            try:
                return DEFAULTS[klass.attrname]
            except KeyError:
                try:
                    return repr(klass._gen_initial_values(1, 1, 1)[0])
                except NotImplementedError:
                    return "No default values"

        def _level(klass: TopologyAttr) -> str:
            if issubclass(klass, AtomAttr):
                level = "atom"
            elif issubclass(klass, ResidueAttr):
                level = "residue"
            elif issubclass(klass, SegmentAttr):
                level = "segment"
            else:
                raise ValueError
            return level

        def _type(klass: TopologyAttr) -> str:
            return klass.dtype  # type: ignore

        self.table_writer = TableWriter(
            filename="generated/topology/defaults.txt",
            sort=True,
            input_items=TOPOLOGY_CLS,
            column_spec=[
                ("Atom", _atom),
                ("AtomGroup", _atomgroup),
                ("default", _default),
                ("level", _level),
                ("type", _type),
            ],
            lines=[],
        )
        self.table_writer.generate_lines_and_write_table()


if __name__ == "__main__":
    TopologyDefaults()
