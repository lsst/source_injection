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

__all__ = ["VisitInjectConnections", "VisitInjectConfig", "VisitInjectTask"]

from typing import cast

from lsst.pex.config import Field
from lsst.pipe.base.connectionTypes import Input, Output

from .inject_base import BaseInjectConfig, BaseInjectConnections, BaseInjectTask


class VisitInjectConnections(  # type: ignore [call-arg]
    BaseInjectConnections,
    dimensions=("instrument", "visit", "detector"),
):
    """Visit-level connections for source injection tasks."""

    visit_summary = Input(
        doc="A visit summary table containing PSF, PhotoCalib and WCS information.",
        name="visit_summary",
        storageClass="ExposureCatalog",
        dimensions=("visit",),
        deferLoad=True,
    )
    input_exposure = Input(
        doc="Exposure to inject synthetic sources into.",
        name="preliminary_visit_image",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    output_exposure = Output(
        doc="Injected Exposure.",
        name="{injected_prefix}preliminary_visit_image",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    output_catalog = Output(
        doc="Catalog of injected sources.",
        name="{injected_prefix}preliminary_visit_image_catalog",
        storageClass="ArrowAstropy",
        dimensions=("instrument", "visit", "detector"),
    )

    def __init__(self, *, config=None):
        config = cast(VisitInjectConfig, config)

        super().__init__(config=config)
        if not config.external_psf and not config.external_photo_calib and not config.external_wcs:
            self.inputs.remove("visit_summary")


class VisitInjectConfig(  # type: ignore [call-arg]
    BaseInjectConfig,
    pipelineConnections=VisitInjectConnections,
):
    """Visit-level configuration for source injection tasks."""

    # Calibrated data options.
    external_psf = Field[bool](
        doc="If True, use the PSF model from a visit summary table. "
        "If False (default), use the PSF model attached to the input exposure.",
        dtype=bool,
        default=False,
    )
    external_photo_calib = Field[bool](
        doc="If True, use the photometric calibration from a visit summary table. "
        "If False (default), use the photometric calibration attached to the input exposure.",
        dtype=bool,
        default=False,
    )
    external_wcs = Field[bool](
        doc="If True, use the astrometric calibration from a visit summary table. "
        "If False (default), use the astrometric calibration attached to the input exposure.",
        dtype=bool,
        default=False,
    )


class VisitInjectTask(BaseInjectTask):
    """Visit-level class for injecting sources into images."""

    _DefaultName = "visitInjectTask"
    ConfigClass = VisitInjectConfig

    def runQuantum(self, butler_quantum_context, input_refs, output_refs):
        inputs = butler_quantum_context.get(input_refs)
        detector_id = inputs["input_exposure"].getDetector().getId()

        try:
            visit_summary = inputs["visit_summary"].get()
        except KeyError:
            # Use internal PSF, PhotoCalib and WCS.
            inputs["psf"] = inputs["input_exposure"].getPsf()
            inputs["photo_calib"] = inputs["input_exposure"].getPhotoCalib()
            inputs["wcs"] = inputs["input_exposure"].getWcs()
        else:
            # Use external PSF, PhotoCalib and WCS.
            detector_summary = visit_summary.find(detector_id)
            if detector_summary:
                inputs["psf"] = detector_summary.getPsf()
                inputs["photo_calib"] = detector_summary.getPhotoCalib()
                inputs["wcs"] = detector_summary.getWcs()
            else:
                raise RuntimeError(f"No record for detector {detector_id} found in visit summary table.")

        input_keys = ["injection_catalogs", "input_exposure", "sky_map", "psf", "photo_calib", "wcs"]
        outputs = self.run(**{key: value for (key, value) in inputs.items() if key in input_keys})
        butler_quantum_context.put(outputs, output_refs)
