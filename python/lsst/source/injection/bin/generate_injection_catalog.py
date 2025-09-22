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

from lsst.daf.butler import Butler

from ..utils import generate_injection_catalog, ingest_injection_catalog
from .source_injection_help_formatter import SourceInjectionHelpFormatter


def build_argparser():
    """Build an argument parser for this script."""
    parser = ArgumentParser(
        description="""Generate a synthetic source injection catalog.

This script generates a synthetic source injection catalog from user supplied
input parameters. The catalog may be printed to screen, ingested into a butler
repository, or written to disk using the astropy Table API.

On-sky source positions are generated using the quasi-random Halton sequence.
By default, the Halton sequence is seeded using the product of the right
ascension and declination limit ranges. This ensures that the same sequence is
always generated for the same limits. This seed may be overridden by the user.

A unique injection ID is generated for each source. The injection ID encodes
two pieces of information: the unique source identification number and the
version number of the source as specified by the ``number`` parameter. To
achieve this, the unique source ID number is multiplied by `10**n` such that
the sum of the multiplied source ID number with the unique repeated version
number will always be unique. For example, an injection catalog with 3 versions
of each source will have injection IDs = 0, 1, 2, 10, 11, 12, 20, 21, 22, etc.
For number = 20, injection IDs = 0, 1, 2, ..., 17, 18, 19, 100, 101, 102, etc.
If number = 1 (the default) then the injection ID will be a simple sequential
list of integers.

An optional butler query for a WCS dataset type may be provided for use in
generating on-sky source positions. If no WCS is provided, the source positions
will be generated using Cartesian geometry.
""",
        formatter_class=SourceInjectionHelpFormatter,
        epilog="More information is available at https://pipelines.lsst.io.",
        add_help=False,
        argument_default=SUPPRESS,
    )
    # General options.
    parser_general = parser.add_argument_group("General Options")
    parser_general.add_argument(
        "-a",
        "--ra-lim",
        type=float,
        help="Right ascension limits of the catalog in degrees.",
        required=True,
        metavar="VALUE",
        nargs=2,
    )
    parser_general.add_argument(
        "-d",
        "--dec-lim",
        type=float,
        help="Declination limits of the catalog in degrees.",
        required=True,
        metavar="VALUE",
        nargs=2,
    )
    parser_general.add_argument(
        "-m",
        "--mag-lim",
        type=float,
        help="The magnitude limits of the catalog in magnitudes.",
        required=False,
        metavar="VALUE",
        nargs=2,
    )
    parser_general.add_argument(
        "-n",
        "--number",
        type=int,
        help="Number of generated parameter combinations. Ignored if density given.",
        metavar="VALUE",
        default=1,
    )
    parser_general.add_argument(
        "-s",
        "--density",
        type=int,
        help="Desired source density (N/deg^2). If given, number option is ignored.",
        metavar="VALUE",
    )
    parser_general.add_argument(
        "-p",
        "--parameter",
        help="An input parameter definition.",
        metavar=("COLNAME VALUE", "VALUE"),
        nargs="+",
        action="append",
    )
    parser_general.add_argument(
        "--seed",
        type=str,
        help="Seed override when generating quasi-random RA/Dec positions.",
        metavar="SEED",
    )

    # Butler options.
    parser_butler = parser.add_argument_group("Butler Options")
    parser_butler.add_argument(
        "-b",
        "--butler-config",
        type=str,
        help="Location of the butler/registry config file.",
        metavar="TEXT",
    )
    # WCS options.
    parser_wcs = parser.add_argument_group("WCS Options")
    parser_wcs.add_argument(
        "-w",
        "--wcs-type-name",
        help="Dataset type containing a `wcs` component for WCS spatial conversions.",
        metavar="TEXT",
    )
    parser_wcs.add_argument(
        "-c",
        "--collections",
        type=str,
        help="Collections to query for dataset types containing a WCS component.",
        metavar="COLL",
    )
    parser_wcs.add_argument(
        "--where",
        type=str,
        help="A string expression similar to an SQL WHERE clause to query for datasets with a WCS component.",
        metavar="COLL",
    )
    # Ingestion options.
    parser_ingest = parser.add_argument_group("Repository Ingestion Options")
    parser_ingest.add_argument(
        "-i",
        "--injection-band",
        type=str,
        help="Band(s) associated with the generated table.",
        metavar=("BAND", "BAND"),
        nargs="+",
    )
    parser_ingest.add_argument(
        "-o",
        "--output-collection",
        type=str,
        help="Name of the output collection to ingest the injection catalog into.",
        metavar="COLL",
    )
    parser_ingest.add_argument(
        "-t",
        "--dataset-type-name",
        type=str,
        help="Output dataset type name for the ingested source injection catalog.",
        metavar="TEXT",
        default="injection_catalog",
    )
    # Options to write table to disk.
    parser_write = parser.add_argument_group("Write to Disk Options")
    parser_write.add_argument(
        "-f",
        "--filename",
        help="Output filename for the generated injection catalog.",
        metavar="TEXT",
    )
    parser_write.add_argument(
        "--format",
        help="Output injection catalog format, overriding automatic format selection.",
        metavar="TEXT",
    )
    parser_write.add_argument(
        "--overwrite",
        help="Overwrite the output file if it already exists.",
        action="store_true",
    )
    # Help.
    parser_misc = parser.add_argument_group("Miscellaneous Options")
    parser_misc.add_argument(
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

    # Validate all butler options are provided.
    butler_config = vars(args).get("butler_config", "")
    wcs_type_name = vars(args).get("wcs_type_name", "")
    collections = vars(args).get("collections", "")
    num_wcs = sum([bool(wcs_type_name), bool(collections)])
    if num_wcs == 1 or (num_wcs == 2 and not butler_config):
        raise RuntimeError("A butler query for WCS requires a butler repo, a dataset type and a collection.")
    injection_band = vars(args).get("injection_band", "")
    output_collection = vars(args).get("output_collection", "")
    dataset_type_name = vars(args).get("dataset_type_name", "")  # Defaults to "injection_catalog".
    num_ingest = sum([bool(injection_band), bool(output_collection)])
    if num_ingest == 1 or (num_ingest == 2 and not butler_config):
        raise RuntimeError("Catalog ingestion requires a butler repo, a band and an output collection.")

    # Parse the input parameters.
    params = {}
    if hasattr(args, "parameter"):
        for param in args.parameter:
            if len(param) < 2:
                raise RuntimeError("Each parameter must be associated with at least one value.")
            name = param[0]
            try:
                values = [float(x) for x in param[1:]]
            except ValueError:
                values = param[1:]
            params[name] = values

    # Get the input WCS.
    if not wcs_type_name:
        wcs = None
        logger.info("No WCS provided, source positions generated using Cartesian geometry.")
    else:
        butler = Butler.from_config(butler_config)
        where = vars(args).get("where", "")  # Optional where query.
        try:
            dataset_ref = list(
                butler.registry.queryDatasets(
                    datasetType=wcs_type_name,
                    collections=collections,
                    where=where,
                    findFirst=True,
                )
            )[0]
        except IndexError:
            raise RuntimeError(f"No {wcs_type_name} dataset type found in {args.collections} for: {where}.")
        ddhandle = butler.getDeferred(wcs_type_name, dataId=dataset_ref.dataId, collections=collections)
        try:
            wcs = ddhandle.get(component="wcs")
        except KeyError:
            raise RuntimeError(f"No WCS component found for {wcs_type_name} dataset type.")
        logger.info("Using WCS in %s for %s.", wcs_type_name, dataset_ref.dataId)

    # Generate the source injection catalog.
    mag_lim = vars(args).get("mag_lim", None)
    density = vars(args).get("density", None)
    seed = vars(args).get("seed", None)
    table = generate_injection_catalog(
        ra_lim=args.ra_lim,
        dec_lim=args.dec_lim,
        mag_lim=mag_lim,
        wcs=wcs,
        number=args.number,
        density=density,
        seed=seed,
        **params,
    )

    # Save table to disk.
    filename = vars(args).get("filename", False)
    file_format = vars(args).get("format", None)
    overwrite = vars(args).get("overwrite", False)
    if filename:
        table.write(filename, format=file_format, overwrite=overwrite)
        logger.info("Written injection catalog to '%s'.", filename)

    # Ingest table into a butler repo.
    if injection_band:
        writeable_butler = Butler.from_config(butler_config, writeable=True)
        for band in injection_band:
            _ = ingest_injection_catalog(
                writeable_butler=writeable_butler,
                table=table,
                band=band,
                output_collection=output_collection,
                dataset_type_name=dataset_type_name,
            )

    # Print the table to stdout.
    if not filename and not injection_band:
        print(table)
