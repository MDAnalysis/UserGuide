"""
Generate groupmethods.txt:

A table of transplanted methods.
"""

from collections.abc import Callable
from typing import Any

import base
from base import TableWriter
from core import TOPOLOGY_CLS
from MDAnalysis.core.groups import GroupBase
from MDAnalysis.core.topologyattrs import TopologyAttr


class TransplantedMethods:
    def __init__(self) -> None:
        def _generate_input_items() -> list[tuple[TopologyAttr, Any]]:
            items = []
            for klass in TOPOLOGY_CLS:
                for name, method in klass.transplants[GroupBase]:
                    items.append((name, klass, method))
            return [x[1:] for x in sorted(items)]

        def _method(klass: TopologyAttr, method: Callable[[TopologyAttr], str]) -> str:
            return base.sphinx_method(method=method)

        def _description(
            klass: TopologyAttr, method: Callable[[TopologyAttr], str]
        ) -> str:
            assert method.__doc__
            return " ".join(method.__doc__.split(".\n")[0].split())

        def _requires(
            klass: TopologyAttr, method: Callable[[TopologyAttr], str]
        ) -> str:
            return klass.attrname  # type: ignore

        input_items = _generate_input_items()
        self.table_writer = TableWriter(
            filename="generated/topology/groupmethods.txt",
            input_items=input_items,
            column_spec=[
                ("Method", _method),
                ("Description", _description),
                ("Requires", _requires),
            ],
            lines=[],
        )
        self.table_writer.generate_lines_and_write_table()


if __name__ == "__main__":
    TransplantedMethods()
