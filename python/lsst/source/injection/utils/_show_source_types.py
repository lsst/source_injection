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

__all__ = ["show_source_types"]

import inspect
import textwrap

import galsim

from ..inject_base import _ALLOWED_SOURCE_TYPES


def show_source_types(
    allowed_source_types: list[str] = _ALLOWED_SOURCE_TYPES,
    wrap_width: int | None = None,
) -> None:
    """Print the signature of each allowed source type.

    Parameters
    ----------
    allowed_source_types : list[str]
        List of allowed source type names.
    wrap_width : int, optional
        Width to wrap the signature text.

    Notes
    -----
    This function prints the signature of each allowed source type,
    wrapping the text to the specified width.
    """
    if wrap_width is None:
        import shutil

        wrap_width = shutil.get_terminal_size().columns
    wrap_width = int(wrap_width)
    for name in sorted(allowed_source_types):
        if name == "Star":
            cls = getattr(galsim, "DeltaFunction")
        elif name == "Trail":
            cls = type("_TempTrail", (), {"__init__": lambda self, trail_length, flux=None: None})
        elif name == "Stamp":
            cls = type("_TempStamp", (), {"__init__": lambda self, stamp, flux=None: None})
        else:
            cls = getattr(galsim, name)
        try:
            params = []
            for p in inspect.signature(cls).parameters.values():
                if p.name == "gsparams":
                    continue
                elif p.name == "flux":
                    params.append("mag=None")
                else:
                    default = f"={p.default}" if p.default is not inspect.Parameter.empty else ""
                    params.append(f"{p.name}{default}")

            sig_text = ", ".join(params)
            wrapped_sig = textwrap.fill(
                f"({sig_text})", width=wrap_width, initial_indent="  ", subsequent_indent="   "
            )
        except (TypeError, ValueError):
            wrapped_sig = "  <no signature available>"

        print(name + ":")
        print(wrapped_sig)
