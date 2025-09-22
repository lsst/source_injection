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

from argparse import SUPPRESS, ArgumentParser

from ..utils import show_source_types
from .source_injection_help_formatter import SourceInjectionHelpFormatter


def build_argparser():
    """Build an argument parser for this script."""
    parser = ArgumentParser(
        description="""Show available source injection types.

This shows all available source types that can be used with the source
injection package, along with their constructor signatures.
""",
        formatter_class=SourceInjectionHelpFormatter,
        epilog="More information is available at https://pipelines.lsst.io.",
        add_help=False,
        argument_default=SUPPRESS,
    )
    parser.add_argument(
        "-w",
        "--wrap-width",
        help="Width to wrap the signature text.",
        default=None,
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
    args = build_argparser().parse_args()
    show_source_types(wrap_width=args.wrap_width)
