#!/usr/bin/env python
"""
Generates:
    - ../formats/format_overview.txt : table of format overview
    - ../formats/coordinate_readers.txt: Coordinate reader information
    - ../formats/classes/*.txt : small tables of parser/reader/writer classes for each format
"""

from collections import defaultdict
from collections.abc import Iterable
from typing import Any, Final, Literal, Type

import base
from base import TableWriter
from core import DESCRIPTIONS
from MDAnalysis import _CONVERTERS, _PARSERS, _READERS, _SINGLEFRAME_WRITERS

HANDLER_T = Literal[
    "Coordinate reader", "Coordinate writer", "Topology parser", "Converter"
]

PAIR_T = tuple[HANDLER_T, dict[str, Any]]

handlers: Iterable[PAIR_T] = [
    ("Coordinate reader", _READERS),
    ("Coordinate writer", _SINGLEFRAME_WRITERS),
    ("Topology parser", _PARSERS),
    ("Converter", _CONVERTERS),
]

FILE_TYPES: dict[str, dict[HANDLER_T, Any]] = defaultdict(dict)

for handler, dct in handlers:
    for format, klass in dct.items():
        if format in {"CHAIN", "MEMORY", "MINIMAL", "NULL"}:
            continue  # get their own pages
        FILE_TYPES[format][handler] = klass


SORTED_FILE_TYPES: Final = sorted(FILE_TYPES.items())

SUCCESS = "\u2713"  # checkmark
FAIL = ""


def _create_key(fmt: str, handlers: dict[HANDLER_T, Type[Any]]) -> str:
    if fmt in DESCRIPTIONS:
        key = fmt
    else:
        key = list(handlers.values())[0].format[0]

        # raise an informative error
        if key not in DESCRIPTIONS:
            key = fmt
    return key


def _file_type(fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str) -> str:
    return base.sphinx_ref(txt=fmt, label=key, suffix="-format")


def _description(fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str) -> str:
    return DESCRIPTIONS[key]


class FormatOverview:
    def __init__(self) -> None:
        def _topology(fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str) -> str:
            if "Topology parser" in handlers:
                return SUCCESS
            return FAIL

        def _coordinates(
            fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str
        ) -> str:
            if "Coordinate reader" in handlers:
                return SUCCESS
            return FAIL

        def _read(fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str) -> str:
            return SUCCESS

        def _write(fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str) -> str:
            if "Coordinate writer" in handlers:
                return SUCCESS
            if "Converter" in handlers:
                return SUCCESS
            return FAIL

        input_items = [
            (format, handlers, _create_key(format, handlers))
            for format, handlers in SORTED_FILE_TYPES
        ]

        self.table_writer = TableWriter(
            filename="formats/format_overview.txt",
            include_table="Table of all supported formats in MDAnalysis",
            column_spec=[
                ("File type", _file_type),
                ("Description", _description),
                ("Topology", _topology),
                ("Coordinates", _coordinates),
                ("Read", _read),
                ("Write", _write),
            ],
            input_items=input_items,
            lines=[],
        )
        self.table_writer.generate_lines_and_write_table()
        self.table_writer.fields["keys"] = list(zip(*input_items))[2]


class CoordinateReaders:
    def __init__(self) -> None:
        def _velocities(
            fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str
        ) -> str:
            if handlers["Coordinate reader"].units.get("velocity"):
                return SUCCESS
            return FAIL

        def _forces(fmt: str, handlers: dict[HANDLER_T, Type[Any]], key: str) -> str:
            if handlers["Coordinate reader"].units.get("force"):
                return SUCCESS
            return FAIL

        input_items = [
            (format, handlers, _create_key(format, handlers))
            for format, handlers in SORTED_FILE_TYPES
            if "Coordinate reader" in handlers
        ]
        self.table_writer = TableWriter(
            filename="formats/coordinate_readers.txt",
            include_table="Table of supported coordinate readers and the information read",
            input_items=input_items,
            column_spec=[
                ("File type", _file_type),
                ("Description", _description),
                ("Velocities", _velocities),
                ("Forces", _forces),
            ],
            lines=[],
        )
        self.table_writer.generate_lines_and_write_table()


class SphinxClasses:
    def __init__(self, fmt: str):
        def _custom_get_lines() -> list[list[str]]:
            lines = []
            for label, klass in sorted(FILE_TYPES[fmt].items()):
                lines.append(
                    [
                        f"**{label}**",
                        base.sphinx_class(klass=klass, tilde=False),
                    ]
                )
            return lines

        lines = _custom_get_lines()
        self.table_writer = TableWriter(
            filename=f"formats/reference/classes/{fmt}.txt",
            lines=lines,
            column_spec=[],
        )
        self.table_writer.write_table()


def main() -> None:
    overview = FormatOverview()
    CoordinateReaders()
    for key in set(overview.table_writer.fields["keys"]):
        SphinxClasses(key)


if __name__ == "__main__":
    main()
