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

import numpy as np
from lsst.pex.config import Field
from lsst.pipe.base.connectionTypes import Input, Output
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.metrics import mean_squared_error

from .inject_base import BaseInjectConfig, BaseInjectConnections, BaseInjectTask


class VisitInjectConnections(  # type: ignore [call-arg]
    BaseInjectConnections,
    dimensions=("instrument", "visit", "detector"),
):
    """Visit-level connections for source injection tasks."""

    visit_summary = Input(
        doc="A visit summary table containing PSF, PhotoCalib and WCS information.",
        name="finalVisitSummary",
        storageClass="ExposureCatalog",
        dimensions=("visit",),
        deferLoad=True,
    )
    input_exposure = Input(
        doc="Exposure to inject synthetic sources into.",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    output_exposure = Output(
        doc="Injected Exposure.",
        name="{injected_prefix}calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
    )
    output_catalog = Output(
        doc="Catalog of injected sources.",
        name="{injected_prefix}calexp_catalog",
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
    brightness_cutoff = Field[float](
        doc="Ignore pixels with image flux below this level when fitting flux vs. variance.",
        dtype=float,
        default=500,
    )
    variance_fit_seed = Field[int](doc="Seed for RANSAC fit of flux vs. variance.", default=0)


class VisitInjectTask(BaseInjectTask):
    """Visit-level class for injecting sources into images."""

    _DefaultName = "visitInjectTask"
    ConfigClass = VisitInjectConfig

    def run(self, injection_catalogs, input_exposure, psf, photo_calib, wcs):
        self.log.info("Fitting flux vs. variance in each pixel.")
        variance_scale = self.get_variance_scale(input_exposure)
        self.log.info("Variance scale factor: %.3f", variance_scale)

        return super().run(injection_catalogs, input_exposure, psf, photo_calib, wcs, variance_scale)

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

    """
    Establish the variance scale by a linear fit of variance vs. flux.
    In practice we see that most pixels in a coadd obey a consistent, simple
    linear relationship between variance and flux, but a small sample skews
    somewhat away from linearity and results in a fit that seems less than
    ideal.

    To identify and hence ignore these odd pixels, we perform a RANSAC fit.
    RANSAC iteratively finds a least-squares straight line fit on random
    subsamples of points, while identifying outliers. The final results tend to
    be much more stable and outlier-robust than a simple linear fit.

    Pixels below a flux level specified by brightness_cutoff are more likely to
    have wild outliers, so for simplicity's sake we simply ignore them.
    """

    def get_variance_scale(self, exposure):
        flux = exposure.image.array.ravel()
        var = exposure.variance.array.ravel()

        # Ignore pixels with nan or infinite values
        good_pixels = np.isfinite(flux) & np.isfinite(var)
        flux = flux[good_pixels]
        var = var[good_pixels]

        # Only fit pixels with at least this much inst flux
        bright_pixels = flux > self.config.brightness_cutoff

        # Simple linear regression to establish MSE
        linear = LinearRegression()
        linear.fit(flux[bright_pixels].reshape(-1, 1), var[bright_pixels])
        linear_mse = mean_squared_error(var, linear.predict(flux.reshape(-1, 1)))

        # RANSAC regression
        ransac = RANSACRegressor(
            loss="squared_error",
            residual_threshold=0.1 * linear_mse,
            random_state=self.config.variance_fit_seed,
        )
        ransac.fit(flux[bright_pixels].reshape(-1, 1), var[bright_pixels])
        variance_scale = float(ransac.estimator_.coef_[0])

        if variance_scale < 0:
            self.log.warning("Slope of final variance vs. flux fit is negative.")

        return variance_scale
