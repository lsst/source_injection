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

__all__ = ["ExposureInjectConnections", "ExposureInjectConfig", "ExposureInjectTask"]

from lsst.pipe.base.connectionTypes import Input, Output

from .inject_visit import VisitInjectConfig, VisitInjectConnections, VisitInjectTask


class ExposureInjectConnections(  # type: ignore [call-arg]
    VisitInjectConnections,
    dimensions=("instrument", "exposure", "detector"),
):
    """Exposure-level connections for source injection tasks."""

    input_exposure = Input(
        doc="Exposure to inject synthetic sources into.",
        name="post_isr_image",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
    )
    output_exposure = Output(
        doc="Injected Exposure.",
        name="{injected_prefix}post_isr_image",
        storageClass="Exposure",
        dimensions=("instrument", "exposure", "detector"),
    )
    output_catalog = Output(
        doc="Catalog of injected sources.",
        name="{injected_prefix}post_isr_image_catalog",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "exposure", "detector"),
    )


class ExposureInjectConfig(  # type: ignore [call-arg]
    VisitInjectConfig,
    pipelineConnections=ExposureInjectConnections,
):
    """Exposure-level configuration for source injection tasks."""

    pass


class ExposureInjectTask(VisitInjectTask):
    """Exposure-level class for injecting sources into images."""

    _DefaultName = "exposureInjectTask"
    ConfigClass = ExposureInjectConfig

    pass
