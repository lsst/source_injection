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

__all__ = ["CoaddInjectConnections", "CoaddInjectConfig", "CoaddInjectTask"]

from lsst.pipe.base.connectionTypes import Input, Output

from .inject_base import BaseInjectConfig, BaseInjectConnections, BaseInjectTask


class CoaddInjectConnections(
    BaseInjectConnections,
    dimensions=("skymap", "tract", "patch", "band"),
    defaultTemplates={
        "coadd_name": "deep",
    },
):
    """Coadd-level connections for source injection tasks."""

    input_exposure = Input(
        doc="Exposure to inject synthetic sources into.",
        name="{coadd_name}Coadd",
        storageClass="ExposureF",
        dimensions=("skymap", "tract", "patch", "band"),
    )
    output_exposure = Output(
        doc="Injected Exposure.",
        name="{injected_prefix}{coadd_name}Coadd",
        storageClass="ExposureF",
        dimensions=("skymap", "tract", "patch", "band"),
    )
    output_catalog = Output(
        doc="Catalog of injected sources.",
        name="{injected_prefix}{coadd_name}Coadd_catalog",
        storageClass="ArrowAstropy",
        dimensions=("skymap", "tract", "patch", "band"),
    )


class CoaddInjectConfig(  # type: ignore [call-arg]
    BaseInjectConfig,
    pipelineConnections=CoaddInjectConnections,
):
    """Coadd-level configuration for source injection tasks."""

    pass


class CoaddInjectTask(BaseInjectTask):
    """Coadd-level class for injecting sources into images."""

    _DefaultName = "coaddInjectTask"
    ConfigClass = CoaddInjectConfig

    def runQuantum(self, butler_quantum_context, input_refs, output_refs):
        inputs = butler_quantum_context.get(input_refs)

        inputs["psf"] = inputs["input_exposure"].getPsf()
        inputs["photo_calib"] = inputs["input_exposure"].getPhotoCalib()
        inputs["wcs"] = inputs["input_exposure"].getWcs()

        input_keys = ["injection_catalogs", "input_exposure", "sky_map", "psf", "photo_calib", "wcs"]
        outputs = self.run(**{key: value for (key, value) in inputs.items() if key in input_keys})
        butler_quantum_context.put(outputs, output_refs)
