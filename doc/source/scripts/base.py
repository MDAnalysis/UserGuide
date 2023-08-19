import collections
import os
import pathlib
import sys
import textwrap
from collections.abc import Callable, Iterable
from typing import Any, Optional, Type

import pandas as pd
import tabulate


def _run_method(method: Callable[..., str], *args: Any) -> str:
    val = method(*args)
    return val


def _generate_row(
    *, column_spec: list[tuple[str, Callable[..., str]]], args: Iterable[Any]
) -> dict[str, str]:
    row = {}
    for heading, method in column_spec:
        val = _run_method(method, *args)
        row[heading] = val
    return row


def _generate_table(
    *,
    input_items: Iterable[Any],
    column_spec: list[tuple[str, Callable[..., str]]],
) -> pd.DataFrame:
    rows = []
    for args in input_items:
        if not isinstance(args, Iterable):
            args = [args]
        line = _generate_row(column_spec=column_spec, args=args)
        rows.append(line)
    df = pd.DataFrame(rows)
    return df


def write_table(
    *,
    path: str,
    headers: list[str],
    lines: list[list[str]],
    include_table: Optional[str] = None,
) -> None:
    parent_directory = pathlib.Path(path).parent
    parent_directory.mkdir(exist_ok=True, parents=True)
    with open(path, "w") as f:
        f.write(f'..\n    Generated by {sys.argv[0].split("/")[-1]}\n\n')
        if include_table:
            f.write(f".. table:: {include_table}\n\n")
        tabled = tabulate.tabulate(lines, headers=headers, tablefmt="rst")
        if include_table:
            tabled = textwrap.indent(tabled, "    ")
        f.write(tabled)
    print("Wrote ", path)


class TableWriter:
    """
    For writing tables with easy column switching.

    Define _set_up_input() and column methods.

    Filename relative to source.
    """

    def __init__(
        self,
        column_spec: list[tuple[str, Callable[..., str]]],
        lines: list[list[str]],
        filename: str = "",
        include_table: Optional[str] = None,
        sort: bool = True,
        input_items: Optional[Iterable[Any]] = None,
    ):
        if column_spec:
            assert input_items
        assert (column_spec and not lines) or (lines and not column_spec)
        stem = os.getcwd().split("source")[0]
        self.path = os.path.join(stem, "source", filename)
        self.filename = filename
        self.include_table = include_table
        self.sort = sort
        self.input_items = input_items
        self.column_spec = column_spec
        self.lines = lines
        self._df = pd.DataFrame()

    @property
    def headers(self) -> list[str]:
        return [column_name for column_name, _ in self.column_spec]

    @property
    def fields(self) -> pd.DataFrame:
        return self._df

    def generate_lines_and_write_table(self) -> None:
        df = _generate_table(
            input_items=self.input_items or [],
            column_spec=self.column_spec,
        )

        lines = df.values.tolist()
        lines = sorted(lines) if self.sort else lines
        self.lines = lines
        self._df = df
        self.write_table()

    def write_table(self) -> None:
        write_table(
            path=self.path,
            headers=self.headers,
            lines=self.lines,
            include_table=self.include_table,
        )


# ==== HELPER FUNCTIONS ==== #


def sphinx_class(*, klass: Type[Any], tilde: bool = True) -> str:
    prefix = "~" if tilde else ""
    return f":class:`{prefix}{klass.__module__}.{klass.__name__}`"


def sphinx_method(*, method: Callable[..., Any], tilde: bool = True) -> str:
    prefix = "~" if tilde else ""
    return ":meth:`{}{}.{}`".format(prefix, method.__module__, method.__qualname__)


def sphinx_ref(*, txt: str, label: Optional[str] = None, suffix: str = "") -> str:
    return f":ref:`{txt} <{label}{suffix}>`"


def sphinx_link(*, txt: str) -> str:
    return "`{}`_".format(txt)
