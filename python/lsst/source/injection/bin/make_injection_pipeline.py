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
import os
import time
from argparse import SUPPRESS, ArgumentParser

from ..utils import make_injection_pipeline
from .source_injection_help_formatter import SourceInjectionHelpFormatter


def build_argparser():
    """Build an argument parser for this script."""
    parser = ArgumentParser(
        description="""Make an expanded source injection pipeline.

This script takes a reference pipeline definition file in YAML format and
prefixes all post-injection dataset type names with the injected prefix. If an
optional injection pipeline definition YAML file is also provided, the
injection task will be merged into the pipeline.

Unless explicitly excluded, all subsets from the reference pipeline containing
the task which generates the injection dataset type will also be updated to
include the injection task. A series of new injection subsets will also be
constructed. These new subsets are copies of existent subsets, but with tasks
not directly impacted by source injection removed. Injected subsets will be the
original existent subset name with the 'injected_' prefix prepended.

When the injection pipeline is constructed, a check on all existing pipeline
contracts is performed. If any contracts are violated, they are removed from
the pipeline. A warning is logged for each contract that is removed.
""",
        formatter_class=SourceInjectionHelpFormatter,
        epilog="More information is available at https://pipelines.lsst.io.",
        add_help=False,
        argument_default=SUPPRESS,
    )
    parser.add_argument(
        "-t",
        "--dataset-type-name",
        type=str,
        help="Name of the datset type being injected into.",
        required=True,
        metavar="TEXT",
    )
    parser.add_argument(
        "-r",
        "--reference-pipeline",
        type=str,
        help="Location of a reference pipeline definition YAML file.",
        required=True,
        metavar="FILE",
    )
    parser.add_argument(
        "-i",
        "--injection-pipeline",
        type=str,
        help="Location of an injection pipeline definition YAML file stub. If "
        "this is not explicitly provided, an attempt to infer the injection "
        "pipeline stub will be made using the injected dataset type name.",
        metavar="FILE",
    )
    parser.add_argument(
        "-e",
        "--exclude-subsets",
        help="Do not update pipeline subsets to include the injection task.",
        action="store_true",
    )
    parser.add_argument(
        "-x",
        "--excluded-tasks",
        type=str,
        help="Comma-separated set of task labels to exclude from the pipeline.",
        metavar="task",
        default="jointcal,gbdesAstrometricFit,fgcmBuildFromIsolatedStars,fgcmFitCycle,fgcmOutputProducts",
    )
    parser.add_argument(
        "-f",
        "--filename",
        help="Path to save a modified pipeline definition YAML file.",
        metavar="FILE",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite the output saved pipeline definition file if it already exists.",
        action="store_true",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="Prefix to prepend to each affected post-injection dataset type name.",
        default="injected_",
    )
    parser.add_argument(
        "--instrument",
        type=str,
        help="Add instrument overrides. Must be a fully qualified class name.",
        metavar="instrument",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Config override for a task, in the format 'label:key=value'.",
        action="append",
    )
    parser.add_argument(
        "-a",
        "--additional-pipelines",
        type=str,
        help="Location(s) of additional input pipeline definition YAML file(s)."
        "Tasks from these additional pipelines will be added to the output injection pipeline.",
        metavar="FILE",
        nargs="+",
    )
    parser.add_argument(
        "-s",
        "--subset-names",
        type=str,
        help="All tasks from any additional pipelines will be added to this subset."
        "The subset will be created if it does not already exist.",
        metavar="FILE",
    )
    parser.add_argument(
        "-d",
        "--subset-description",
        type=str,
        help="The description given to a new subset which holds tasks from additional pipelines provided."
        "Note: this argument is ignored if the subset already exists.",
        metavar="FILE",
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
        format="%(levelname)s %(asctime)s.%(msecs)03d" + tz + " - %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    args = build_argparser().parse_args()
    if hasattr(args, "filename"):
        if os.path.exists(args.filename):
            if not hasattr(args, "overwrite"):
                raise RuntimeError(f"File {args.filename} already exists; use --overwrite to write anyway.")
            else:
                logger.warning("File %s already exists; overwriting.", args.filename)
        pipeline = make_injection_pipeline(
            **{k: v for k, v in vars(args).items() if k not in ["filename", "overwrite"]}
        )
        pipeline.write_to_uri(args.filename)
        logger.info(
            "Modified pipeline definition YAML file saved at %s.",
            os.path.realpath(args.filename),
        )
    else:
        pipeline = make_injection_pipeline(
            **{k: v for k, v in vars(args).items() if k not in ["filename", "overwrite"]}
        )
        print("\n", pipeline, sep="")
