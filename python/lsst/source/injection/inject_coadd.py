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
    max_trials_1 = Field[int](
        doc="Maximum number of trials the first RANSAC fit is allowed to run.",
        default=1000,
    )
    max_trials_2 = Field[int](
        doc="Maximum number of trials the second RANSAC fit is allowed to run.",
        default=1000,
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
    kmeans_n_init_1 = Field[int](
        doc="Number of times the first k-means clustering is run with different initial centroids.",
        default=10,
    )
    kmeans_n_init_2 = Field[int](
        doc="Number of times the second k-means clustering is run with different initial centroids.",
        default=10,
    )


class CoaddInjectTask(BaseInjectTask):
    """Coadd-level class for injecting sources into images."""

    _DefaultName = "coaddInjectTask"
    ConfigClass = CoaddInjectConfig

    def run(self, injection_catalogs, input_exposure, psf, photo_calib, wcs):
        self.log.info("Fitting flux vs. variance in each pixel.")
        variance_scale = self.get_variance_scale(input_exposure)
        self.log.info("Variance scale factor: %.6f", variance_scale)

        return super().run(injection_catalogs, input_exposure, psf, photo_calib, wcs, variance_scale)

    def runQuantum(self, butler_quantum_context, input_refs, output_refs):
        inputs = butler_quantum_context.get(input_refs)

        inputs["psf"] = inputs["input_exposure"].getPsf()
        inputs["photo_calib"] = inputs["input_exposure"].getPhotoCalib()
        inputs["wcs"] = inputs["input_exposure"].getWcs()

        input_keys = ["injection_catalogs", "input_exposure", "sky_map", "psf", "photo_calib", "wcs"]
        outputs = self.run(**{key: value for (key, value) in inputs.items() if key in input_keys})
        butler_quantum_context.put(outputs, output_refs)

    """
    Establish the variance scale by a linear fit of variance vs. flux.
    In practice we see that most pixels in a coadd obey a consistent, simple
    linear relationship between variance and flux, but a small sample do not.

    To identify and hence ignore these odd pixels, we perform two rounds of
    RANSAC fits. RANSAC iteratively finds a least-squares straight line fit on
    random subsamples of points, while identifying outliers. The inlier pixels
    from a first RANSAC fit tend to be regions of empty space and the extreme
    outer edges of galaxies, while the outliers tend to be the inner regions of
    galaxies.

    We can pick up some more pixels in the galaxies by running a second RANSCAC
    fit on just the outliers of the first fit. In the second round, the inliers
    tend to be the bulk of the galaxy, while the outliers are the innermost
    cores of galaxies. In practice, the inliers from this second round tend to
    have a qualitatively similar variance-vs-flux relationship to the outliers
    from the first fit and do not strongly alter the final variance scale we
    get; while the outliers from the second round are truly wild and should
    clearly be omitted from a simple linear fit.

    In both rounds, random variation of the points sampled by the RANSAC fit
    causes the fitted slope and intercept to vary, and sometimes the fit can
    settle on a pathological sample of points as its inliers. We try to avoid
    such pathologies by running each fit multiple times with different seeds,
    and using K-Means clustering on the resulting slopes to identify the most
    stable value.
    """

    def get_variance_scale(self, exposure):
        flux = exposure.image.array.ravel()
        var = exposure.variance.array.ravel()

        # Ignore pixels with nan or infinite values
        good_pixels = np.isfinite(flux) & np.isfinite(var)
        flux = flux[good_pixels]
        var = var[good_pixels]

        # Simple linear regression to establish MSE.
        linear = LinearRegression()
        linear.fit(flux.reshape(-1, 1), var)
        linear_mse = mean_squared_error(var, linear.predict(flux.reshape(-1, 1)))

        # First RANSAC fit
        fit_results = []
        for seed in range(
            self.config.variance_fit_seed_1, self.config.variance_fit_seed_1 + self.config.n_fits_1
        ):
            ransac = RANSACRegressor(
                loss="squared_error",
                residual_threshold=self.config.threshold_scale_1 * linear_mse,
                max_trials=self.config.max_trials_1,
                random_state=seed,
            )
            ransac.fit(flux.reshape(-1, 1), var)
            # Remember fit results
            slope = ransac.estimator_.coef_[0]
            fit_results.append((slope, seed))

        # K-means cluster the first round of fits,
        # to find the most stable results.
        kmeans = KMeans(
            n_clusters=self.config.n_clusters_1,
            random_state=self.config.kmeans_seed_1,
            n_init=self.config.kmeans_n_init_1,
        )
        kmeans.fit(np.log(np.array([f[0] for f in fit_results if f[0] > 0])).reshape(-1, 1))
        label_counts = [np.sum(kmeans.labels_ == idx) for idx in range(self.config.n_clusters_1)]

        # Recall one of the fits, chosen arbitrarily from those which are both
        # stable, according to the first k-means fit, and positive-slope.
        stable_fit_seeds = np.array([f[1] for f in fit_results if f[0] > 0])[
            kmeans.labels_ == np.argmax(label_counts)
        ]
        if len(stable_fit_seeds == 0):
            # No positive-slope fit found.
            # Allow the fitted slope to be negative but throw a warning.
            self.log.warning(
                "No positive-slope result in the first round of "
                "RANSAC fits. Proceeding with a negative-slope fit."
            )
            stable_fit_seeds = np.array([f[1] for f in fit_results])[
                kmeans.labels_ == np.argmax(label_counts)
            ]
        seed = stable_fit_seeds[0]
        ransac = RANSACRegressor(
            loss="squared_error",
            residual_threshold=self.config.threshold_scale_1 * linear_mse,
            max_trials=self.config.max_trials_1,
            random_state=seed,
        )
        ransac.fit(flux.reshape(-1, 1), var)

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
                max_trials=self.config.max_trials_2,
                random_state=seed,
            )
            ransac.fit(flux[outlier_mask_1].reshape(-1, 1), var[outlier_mask_1])
            # Remember fit results
            slope = ransac.estimator_.coef_[0]
            fit_results.append((slope, seed))

        # K-Means cluster the second round of fits,
        # to find the most stable result
        kmeans = KMeans(
            n_clusters=self.config.n_clusters_2,
            random_state=self.config.kmeans_seed_2,
            n_init=self.config.kmeans_n_init_2,
        )
        kmeans.fit(np.log(np.array([f[0] for f in fit_results if f[0] > 0])).reshape(-1, 1))
        label_counts = [np.sum(kmeans.labels_ == idx) for idx in range(self.config.n_clusters_2)]

        # Recall one of the stable fits
        stable_fit_seeds = np.array([f[1] for f in fit_results if f[0] > 0])[
            kmeans.labels_ == np.argmax(label_counts)
        ]
        if len(stable_fit_seeds == 0):
            # No positive-slope fit found.
            # Allow the fitted slope to be negative but throw a warning.
            self.log.warning(
                "No positive-slope result in the second round of "
                "RANSAC fits. Proceeding with a negative-slope fit."
            )
            stable_fit_seeds = np.array([f[1] for f in fit_results])[
                kmeans.labels_ == np.argmax(label_counts)
            ]
        seed = stable_fit_seeds[0]
        ransac = RANSACRegressor(
            loss="squared_error",
            residual_threshold=self.config.threshold_scale_2 * linear_mse,
            max_trials=self.config.max_trials_2,
            random_state=seed,
        )
        ransac.fit(flux[outlier_mask_1].reshape(-1, 1), var[outlier_mask_1])

        # Pixels with a "good" variance vs. flux relationship:
        # Union of the inliers from the first fit
        # together with the inliers from the second fit.
        flux_good = np.concatenate(
            (flux[inlier_mask_1], flux[outlier_mask_1][ransac.inlier_mask_]),
            axis=None,
        )
        var_good = np.concatenate(
            (var[inlier_mask_1], var[outlier_mask_1][ransac.inlier_mask_]),
            axis=None,
        )

        # Fit all the good pixels with a simple least squares regression.
        linear = LinearRegression()
        linear.fit(flux_good.reshape(-1, 1), var_good)
        variance_scale = float(linear.coef_[0])

        if variance_scale < 0:
            self.log.warning("Slope of final variance vs. flux fit is negative.")

        # Return the slope of the final fit.
        return variance_scale
