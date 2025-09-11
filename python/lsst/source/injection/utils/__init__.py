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

"""
The source_injection package provides functionality for injecting synthetic
sources into astronomical images.

This sub-package contains utility functions and tasks for assisting in
generating synthetic source inputs.

Modules
-------
_make_injection_pipeline : Create source injection pipelines.
_generate_injection_catalog : Generate source injection catalogs.
_ingest_injection_catalog : Ingest source injection catalogs into a repository.
_consolidate_injected_coadd_catalogs : Consolidate per-patch injected
    coadd catalogs into a single per-tract table.
"""

from ._consolidate_injected_coadd_catalogs import *  # noqa: F401,F403
from ._generate_injection_catalog import *  # noqa: F401,F403
from ._ingest_injection_catalog import *  # noqa: F401,F403
from ._make_injection_pipeline import *  # noqa: F401,F403
from ._show_source_types import *  # noqa: F401,F403
