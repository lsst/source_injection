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

import numpy as np
from lsst.pex.config import Field
from lsst.pipe.base.connectionTypes import Input, Output
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.metrics import mean_squared_error

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

    n_fits_1 = Field[int](
        doc="Perform this many RANSAC fits in the first round, to get a sample "
        "of different slopes based on the different random samples of points used "
        "in the fit.",
        default=20,
    )
    n_fits_2 = Field[int](
        doc="Perform this many RANSAC fits in the second round, to get a sample "
        "of different slopes based on the different random samples of points used "
        "in the fit.",
        default=20,
    )
    threshold_scale_1 = Field[float](
        doc="An outlier in the first RANSAC fit is farther from the fit line, "
        "in terms of squared error, than this multiple of the initial linear MSE.",
        default=0.1,
    )
    threshold_scale_2 = Field[float](
        doc="An outlier in the second RANSAC fit is farther from the fit line, "
        "in terms of squared error, than this multiple of the initial linear MSE.",
        default=0.1,
    )
    variance_fit_seed_1 = Field[int](doc="Seed for first RANSAC fit of flux vs. variance.", default=0)
    variance_fit_seed_2 = Field[int](doc="Seed for second RANSAC fit of flux vs. variance.", default=0)
    n_clusters_1 = Field[int](
        doc="K-means cluster the first set of RANSAC fits using this many clusters, "
        "in order to find the most stable slope (biggest cluster).",
        default=4,
    )
    n_clusters_2 = Field[int](
        doc="K-means cluster the second set of RANSAC fits using this many clusters, "
        "in order to find the most stable slope (biggest cluster).",
        default=3,
    )
    kmeans_seed_1 = Field[int](doc="Seed for first round of k-means clustering.", default=0)
    kmeans_seed_2 = Field[int](doc="Seed for second round of k-means clustering.", default=0)


class CoaddInjectTask(BaseInjectTask):
    """Coadd-level class for injecting sources into images."""

    _DefaultName = "coaddInjectTask"
    ConfigClass = CoaddInjectConfig

    def run(self, injection_catalogs, input_exposure, psf, photo_calib, wcs):
        self.log.info("Fitting flux vs. variance in each pixel.")
        self.config.variance_scale = self.get_variance_scale(input_exposure)
        self.log.info("Variance scale factor: %.6f", self.config.variance_scale)

        return super().run(injection_catalogs, input_exposure, psf, photo_calib, wcs)

    def runQuantum(self, butler_quantum_context, input_refs, output_refs):
        inputs = butler_quantum_context.get(input_refs)

        inputs["psf"] = inputs["input_exposure"].getPsf()
        inputs["photo_calib"] = inputs["input_exposure"].getPhotoCalib()
        inputs["wcs"] = inputs["input_exposure"].getWcs()

        input_keys = ["injection_catalogs", "input_exposure", "sky_map", "psf", "photo_calib", "wcs"]
        outputs = self.run(**{key: value for (key, value) in inputs.items() if key in input_keys})
        butler_quantum_context.put(outputs, output_refs)

    def get_variance_scale(self, exposure):
        x = exposure.image.array.flatten()
        y = exposure.variance.array.flatten()

        # Identify bad pixels
        bad_pixels = ~np.isfinite(x) | ~np.isfinite(y)
        # Replace bad pixel values with the image median
        if np.sum(bad_pixels) > 0:
            median_image_value = np.median(x)
            median_variance_value = np.median(y)
            x[bad_pixels] = median_image_value
            y[bad_pixels] = median_variance_value

        # Simple linear regression to establish MSE.
        linear = LinearRegression()
        linear.fit(x.reshape(-1, 1), y)
        linear_mse = mean_squared_error(y, linear.predict(x.reshape(-1, 1)))

        # First RANSAC fit
        fit_results = []
        for seed in range(
            self.config.variance_fit_seed_1, self.config.variance_fit_seed_1 + self.config.n_fits_1
        ):
            ransac = RANSACRegressor(
                loss="squared_error",
                residual_threshold=self.config.threshold_scale_1 * linear_mse,
                max_trials=1000,
                random_state=seed,
            )
            ransac.fit(x.reshape(-1, 1), y)
            # Remember fit results
            slope = ransac.estimator_.coef_[0]
            fit_results.append((slope, seed))

        # K-means cluster the first round of fits,
        # to find the most stable results.
        kmeans = KMeans(
            n_clusters=self.config.n_clusters_1, random_state=self.config.kmeans_seed_1, n_init=10
        )
        kmeans.fit(np.log(np.array([f[0] for f in fit_results if f[0] > 0])).reshape(-1, 1))
        label_counts = [np.sum(kmeans.labels_ == idx) for idx in range(self.config.n_clusters_1)]

        # Recall one of the fits, chosen arbitrarily from those which are both
        # stable, according to the first k-means fit, and positive-slope.
        stable_fit_seeds = np.array([f[1] for f in fit_results if f[0] > 0])[
            kmeans.labels_ == np.argmax(label_counts)
        ]
        if len(stable_fit_seeds == 0):
            # Throw a warning
            pass
        else:
            seed = stable_fit_seeds[0]
        ransac = RANSACRegressor(
            loss="squared_error",
            residual_threshold=self.config.threshold_scale_1 * linear_mse,
            max_trials=1000,
            random_state=seed,
        )
        ransac.fit(x.reshape(-1, 1), y)

        # Label the pixels with a "good" variance vs. flux relationship
        # (the inliers), together with the ones that are further from a simple
        # straight line (the outliers).
        inlier_mask_1 = ransac.inlier_mask_
        outlier_mask_1 = ~inlier_mask_1

        # Second RANSAC fit,
        # on just the outliers of the 1st fit.
        fit_results = []
        for seed in range(
            self.config.variance_fit_seed_2, self.config.variance_fit_seed_2 + self.config.n_fits_2
        ):
            ransac = RANSACRegressor(
                loss="squared_error",
                residual_threshold=self.config.threshold_scale_2 * linear_mse,
                max_trials=1000,
                random_state=seed,
            )
            ransac.fit(x[outlier_mask_1].reshape(-1, 1), y[outlier_mask_1])
            # Remember fit results
            slope = ransac.estimator_.coef_[0]
            fit_results.append((slope, seed))

        # K-Means cluster the second round of fits,
        # to find the most stable result
        kmeans = KMeans(
            n_clusters=self.config.n_clusters_2, random_state=self.config.kmeans_seed_2, n_init=10
        )
        kmeans.fit(np.log(np.array([f[0] for f in fit_results if f[0] > 0])).reshape(-1, 1))
        label_counts = [np.sum(kmeans.labels_ == idx) for idx in range(self.config.n_clusters_2)]

        # Recall one of the stable fits
        stable_fit_seeds = np.array([f[1] for f in fit_results if f[0] > 0])[
            kmeans.labels_ == np.argmax(label_counts)
        ]
        if len(stable_fit_seeds == 0):
            # Throw a warning
            pass
        else:
            seed = stable_fit_seeds[0]
        ransac = RANSACRegressor(
            loss="squared_error",
            residual_threshold=self.config.threshold_scale_2 * linear_mse,
            max_trials=1000,
            random_state=seed,
        )
        ransac.fit(x[outlier_mask_1].reshape(-1, 1), y[outlier_mask_1])

        # Pixels with a "good" variance vs. flux relationship:
        # Union of the inliers from the first fit
        # together with the inliers from the second fit.
        x_good = np.concatenate((x[inlier_mask_1], x[outlier_mask_1][ransac.inlier_mask_]), axis=None)
        y_good = np.concatenate((y[inlier_mask_1], y[outlier_mask_1][ransac.inlier_mask_]), axis=None)

        # Fit all the good pixels with a simple least squares regression.
        linear = LinearRegression()
        linear.fit(x_good.reshape(-1, 1), y_good)

        # Return the slope of the final fit.
        return float(linear.coef_[0])
