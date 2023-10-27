#!/usr/bin/env python
"""
Generate three tables:
    - connectivityattrs.txt: A table of supported formats for bonds, angles, dihedrals, impropers
    - topologyattrs.txt: A table of supported formats for non-connectivity attributes.
    - topology_parsers.txt: A table of all formats and the attributes they read and guess

This script imports the testsuite, which tests these.
"""
from collections import defaultdict
from typing import Any

import base
from base import TableWriter
from core import DESCRIPTIONS, NON_CORE_ATTRS
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysisTests.topology.base import mandatory_attrs
from MDAnalysisTests.topology.test_crd import TestCRDParser
from MDAnalysisTests.topology.test_dlpoly import (
    TestDLPConfigParser,
    TestDLPHistoryParser,
)
from MDAnalysisTests.topology.test_dms import TestDMSParser
from MDAnalysisTests.topology.test_fhiaims import TestFHIAIMS
from MDAnalysisTests.topology.test_gms import GMSBase
from MDAnalysisTests.topology.test_gro import TestGROParser
from MDAnalysisTests.topology.test_gsd import TestGSDParser
from MDAnalysisTests.topology.test_hoomdxml import TestHoomdXMLParser
from MDAnalysisTests.topology.test_lammpsdata import LammpsBase, TestDumpParser
from MDAnalysisTests.topology.test_mmtf import TestMMTFParser
from MDAnalysisTests.topology.test_mol2 import TestMOL2Base
from MDAnalysisTests.topology.test_pdb import TestPDBParser
from MDAnalysisTests.topology.test_pdbqt import TestPDBQT
from MDAnalysisTests.topology.test_pqr import TestPQRParser
from MDAnalysisTests.topology.test_psf import PSFBase
from MDAnalysisTests.topology.test_top import TestPRMParser
from MDAnalysisTests.topology.test_tprparser import TPRAttrs
from MDAnalysisTests.topology.test_txyz import TestTXYZParser
from MDAnalysisTests.topology.test_xpdb import TestXPDBParser
from MDAnalysisTests.topology.test_xyz import XYZBase

PARSER_TESTS = (
    TestCRDParser,
    TestDLPHistoryParser,
    TestDLPConfigParser,
    TestDMSParser,
    TestFHIAIMS,
    GMSBase,
    TestGROParser,
    TestGSDParser,
    TestHoomdXMLParser,
    LammpsBase,
    TestMMTFParser,
    TestMOL2Base,
    TestPDBParser,
    TestPDBQT,
    TestPQRParser,
    PSFBase,
    TestPRMParser,
    TPRAttrs,
    TestTXYZParser,
    TestXPDBParser,
    XYZBase,
    TestDumpParser,
)


def create_parser_attributes() -> dict[Any, tuple[set[str], set[str]]]:
    parser_attrs = {}
    for test_parser_class in PARSER_TESTS:
        expected = set(test_parser_class.expected_attrs) - set(mandatory_attrs)
        guessed = set(test_parser_class.guessed_attrs)
        # clunky hack for PDB
        if test_parser_class is TestPDBParser:
            expected.add("elements")
        parser_attrs[test_parser_class.parser] = (expected, guessed)
    return parser_attrs


class TopologyParsers:
    def __init__(self) -> None:
        def _keys(parser: TopologyReaderBase) -> tuple[str, str]:
            f = parser.format
            if isinstance(f, (list, tuple)):
                key = f[0]
                label = ", ".join(f)
            else:
                key = label = f
            return (key, label)

        def _description(
            parser: TopologyReaderBase,
            expected: set[str],
            guessed: set[str],
            key_label: tuple[str, str],
        ) -> str:
            key, label = key_label
            return DESCRIPTIONS[key]

        def _format(
            parser: TopologyReaderBase,
            expected: set[str],
            guessed: set[str],
            key_label: tuple[str, str],
        ) -> str:
            key, label = key_label
            return base.sphinx_ref(txt=label, label=key, suffix="-format")

        def _attributes_read(
            parser: TopologyReaderBase,
            expected: set[str],
            guessed: set[str],
            key_label: tuple[str, str],
        ) -> str:
            vals = sorted(expected - guessed)
            return ", ".join(vals)

        def _attributes_guessed(
            parser: TopologyReaderBase,
            expected: set[str],
            guessed: set[str],
            key_label: tuple[str, str],
        ) -> str:
            return ", ".join(sorted(guessed))

        parser_attrs = create_parser_attributes()
        input_items = [
            [parser, expected, guessed, _keys(parser=parser)]
            for parser, (expected, guessed) in parser_attrs.items()
        ]
        self.table_writer = TableWriter(
            filename="formats/topology_parsers.txt",
            include_table="Table of supported topology parsers and the attributes read",
            sort=True,
            input_items=input_items,
            lines=[],
            column_spec=[
                ("Format", _format),
                ("Description", _description),
                ("Attributes read", _attributes_read),
                ("Attributes guessed", _attributes_guessed),
            ],
        )
        self.table_writer.generate_lines_and_write_table()


def get_format_attrs(topology_parsers: TopologyParsers) -> dict[str, set[str]]:
    attrs = defaultdict(set)
    writer = topology_parsers.table_writer
    assert writer.input_items
    for format, (_, expected, guessed, _) in zip(
        writer.fields["Format"],
        writer.input_items,
    ):
        for attribute in expected | guessed:
            attrs[attribute].add(format)
    return attrs


class TopologyAttrs:
    def __init__(self, attrs: dict[str, set[str]]) -> None:
        def _atom(name: str, singular: str, description: str) -> str:
            return singular

        def _atomgroup(name: str, singular: str, description: str) -> str:
            return name

        def _description(name: str, singular: str, description: str) -> str:
            return description

        def _supported_formats(
            name: str, singular: str, description: str
        ) -> str:
            return ", ".join(sorted(attrs[name]))

        input_items = sorted(
            [x, *y]
            for x, y in NON_CORE_ATTRS.items()
            if x not in set(mandatory_attrs)
        )
        self.table_writer = TableWriter(
            filename="generated/topology/topologyattrs.txt",
            lines=[],
            input_items=input_items,
            column_spec=[
                ("Atom", _atom),
                ("AtomGroup", _atomgroup),
                ("Description", _description),
                ("Supported formats", _supported_formats),
            ],
        )
        self.table_writer.generate_lines_and_write_table()


class ConnectivityAttrs:
    def __init__(self, attrs: dict[str, set[str]]) -> None:
        def _atom(name: str) -> str:
            return name

        def _atomgroup(name: str) -> str:
            return name

        def _supported_formats(name: str) -> str:
            return ", ".join(sorted(attrs[name]))

        input_items = [("bonds",), ("angles",), ("dihedrals",), ("impropers",)]

        self.table_writer = TableWriter(
            filename="generated/topology/connectivityattrs.txt",
            input_items=input_items,
            lines=[],
            column_spec=[
                ("Atom", _atom),
                ("AtomGroup", _atomgroup),
                ("Supported formats", _supported_formats),
            ],
        )
        self.table_writer.generate_lines_and_write_table()


def main() -> None:
    top = TopologyParsers()
    topology_attrs = get_format_attrs(top)
    TopologyAttrs(topology_attrs)
    ConnectivityAttrs(topology_attrs)


if __name__ == "__main__":
    main()
