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

__all__ = ["SourceInjectionHelpFormatter"]

from argparse import ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter


class SourceInjectionHelpFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
    """Custom help formatter class to override standard formatting output."""

    def __init__(self, prog):
        super().__init__(prog, width=79, max_help_position=6)

    def _format_action_invocation(self, action):
        """Modify help message format to remove duplicate metavar entries."""
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ", ".join(action.option_strings) + " " + args_string

    def _split_lines(self, text, width):
        """Split text into lines; add a line break after each option block."""
        return super()._split_lines(text, width) + [""]
