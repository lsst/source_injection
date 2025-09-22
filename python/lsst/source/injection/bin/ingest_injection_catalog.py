# This file is part of source_injection.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

import logging
import time
from argparse import SUPPRESS, ArgumentParser

from astropy.table import Table

from lsst.daf.butler import Butler
from lsst.daf.butler.formatters.parquet import ParquetFormatter, arrow_to_astropy, pq

from ..utils import ingest_injection_catalog
from .source_injection_help_formatter import SourceInjectionHelpFormatter


def _is_parquet(filename: str):
    """Return if a filename has a parquet extension.

    Notes
    -----
    This could be replaced with astropy.io.misc.parquet_identify, which has
    additional functionality to open the file and validate the first set of
    bytes.
    """
    extensions = {".parquet", ".parq"} | {
        ParquetFormatter.default_extension,
    }
    extensions = tuple(ext for ext in extensions if ext is not None)

    return filename.endswith(extensions)


def build_argparser():
    """Build an argument parser for this script."""
    parser = ArgumentParser(
        description="""Ingest a source injection catalog into the butler.

This script reads an on-disk input catalog or multiple per-band input catalogs
and ingests these data into the butler. Catalogs are read in using the astropy
Table API. Parquet filenames will be read through daf_butler functions, whereas
other file types will attempt to use astropy.table.Table.read. See DM-44159 for
details.

An attempt at auto-identification of the input catalog file format type will be
made for each input. A manually specified format can instead be specified for
all input catalogs using the ``--format`` option.

Each injection catalog must be associated with at least one band, specified
immediately after the path to the input catalog. Multiple space-separated bands
can be provided for a single input catalog. The injection catalog option may
also be called multiple times to ingest multiple per-band catalogs.
""",
        formatter_class=SourceInjectionHelpFormatter,
        epilog="More information is available at https://pipelines.lsst.io.",
        add_help=False,
        argument_default=SUPPRESS,
    )
    parser.add_argument(
        "-b",
        "--butler-config",
        type=str,
        help="Location of the butler/registry config file.",
        required=True,
        metavar="TEXT",
    )
    parser.add_argument(
        "-i",
        "--injection-catalog",
        type=str,
        help="Location of the input source injection catalog and all associated bands.",
        required=True,
        metavar=("FILE BAND", "BAND"),
        nargs="+",
        action="append",
    )
    parser.add_argument(
        "-o",
        "--output-collection",
        type=str,
        help="Name of the output collection to ingest the injection catalog into.",
        required=True,
        metavar="COLL",
    )
    parser.add_argument(
        "-t",
        "--dataset-type-name",
        type=str,
        help="Output dataset type name for the ingested source injection catalog.",
        metavar="TEXT",
        default="injection_catalog",
    )
    parser.add_argument(
        "--format",
        type=str,
        help="Format of the input injection catalog(s), overriding auto-identification.",
        metavar="TEXT",
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.",
    )
    return parser


def main():
    """Use this as the main entry point when calling from the command line."""
    # Set up logging.
    tz = time.strftime("%z")
    logging.basicConfig(
        format="%(levelname)s %(asctime)s.%(msecs)03d" + tz + " - %(message)s", datefmt="%Y-%m-%dT%H:%M:%S"
    )
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    args = build_argparser().parse_args()

    injection_catalogs = []
    for injection_catalog_schema in args.injection_catalog:
        if len(injection_catalog_schema) < 2:
            raise RuntimeError("Each injection catalog must be associated with at least one band.")
        injection_catalog = injection_catalog_schema.pop(0)
        for band in injection_catalog_schema:
            injection_catalogs.append((injection_catalog, band))

    writeable_butler = Butler.from_config(args.butler_config, writeable=True)
    injection_catalog_format = vars(args).get("format", None)

    injection_catalogs_table = Table(rows=injection_catalogs, names=("injection_catalog", "band"))
    injection_catalogs_groups = injection_catalogs_table.group_by("band")

    for injection_catalogs_group in injection_catalogs_groups.groups:
        band = str(injection_catalogs_group["band"][0])

        injection_catalogs_band = []
        for injection_catalog in injection_catalogs_group["injection_catalog"]:
            # The character_as_bytes=False option is preferred, if possible.
            if isinstance(injection_catalog, str) and _is_parquet(injection_catalog):
                tbl = arrow_to_astropy(pq.read_table(injection_catalog, use_threads=False))
            else:
                try:
                    tbl = Table.read(
                        injection_catalog,
                        format=injection_catalog_format,
                        character_as_bytes=False,
                    )
                except TypeError:
                    tbl = Table.read(injection_catalog, format=injection_catalog_format)
            injection_catalogs_band.append(tbl)

        _ = ingest_injection_catalog(
            writeable_butler=writeable_butler,
            table=injection_catalogs_band,
            band=band,
            **{
                k: v
                for k, v in vars(args).items()
                if k not in ["butler_config", "injection_catalog", "format"]
            },
        )
