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


import numpy as np
from lsst.afw.geom import makeCdMatrix, makeSkyWcs
from lsst.afw.image import makePhotoCalibFromCalibZeroPoint
from lsst.geom import Box2I, Extent2I, Point2D, Point2I, SpherePoint, degrees
from lsst.ip.isr.isrTask import IsrTask
from lsst.meas.algorithms.testUtils import plantSources
from lsst.pipe.base import Pipeline
from lsst.pipe.base.pipelineIR import LabeledSubset
from lsst.pipe.tasks.calibrate import CalibrateTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
from lsst.source.injection import generate_injection_catalog


def make_test_exposure():
    """Make a test exposure with a PSF attached and stars placed randomly.

    This function generates a noisy synthetic image with Gaussian PSFs injected
    into the frame.
    The exposure is returned with a WCS, PhotoCalib and PSF attached.

    Returns
    -------
    exposure : `lsst.afw.image.Exposure`
        Exposure with calibs attached and stars placed randomly.
    """
    # Inspired by meas_algorithms test_dynamicDetection.py.
    xy0 = Point2I(12345, 67890)  # xy0 for image
    dims = Extent2I(2345, 2345)  # Dimensions of image
    bbox = Box2I(xy0, dims)  # Bounding box of image
    sigma = 3.21  # PSF sigma
    buffer = 4.0  # Buffer for star centers around edge
    n_sigma = 5.0  # Number of PSF sigmas for kernel
    sky = 12345.6  # Initial sky level
    num_stars = 100  # Number of stars
    noise = np.sqrt(sky) * np.pi * sigma**2  # Poisson noise per PSF
    faint = 1.0 * noise  # Faintest level for star fluxes
    bright = 100.0 * noise  # Brightest level for star fluxes
    star_bbox = Box2I(bbox)  # Area on image in which we can put star centers
    star_bbox.grow(-int(buffer * sigma))  # Shrink star_bbox
    pixel_scale = 1.0e-5 * degrees  # Pixel scale (1E-5 deg = 0.036 arcsec)

    # Make an exposure with a PSF attached; place stars randomly.
    rng = np.random.default_rng(12345)
    stars = [
        (xx, yy, ff, sigma)
        for xx, yy, ff in zip(
            rng.uniform(star_bbox.getMinX(), star_bbox.getMaxX(), num_stars),
            rng.uniform(star_bbox.getMinY(), star_bbox.getMaxY(), num_stars),
            np.linspace(faint, bright, num_stars),
        )
    ]
    exposure = plantSources(bbox, 2 * int(n_sigma * sigma) + 1, sky, stars, True)

    # Set WCS and PhotoCalib.
    exposure.setWcs(
        makeSkyWcs(
            crpix=Point2D(0, 0),
            crval=SpherePoint(0, 0, degrees),
            cdMatrix=makeCdMatrix(scale=pixel_scale),
        )
    )
    exposure.setPhotoCalib(makePhotoCalibFromCalibZeroPoint(1e10, 1e8))
    return exposure


def make_test_injection_catalog(wcs, bbox):
    """Make a test source injection catalog.

    This function generates a test source injection catalog consisting of 30
    star-like sources of varying magnitude.

    Parameters
    ----------
    wcs : `lsst.afw.geom.SkyWcs`
        WCS associated with the exposure.
    bbox : `lsst.geom.Box2I`
        Bounding box of the exposure.

    Returns
    -------
    injection_catalog : `astropy.table.Table`
        Source injection catalog.
    """
    radec0 = wcs.pixelToSky(bbox.getBeginX(), bbox.getBeginY())
    radec1 = wcs.pixelToSky(bbox.getEndX(), bbox.getEndY())
    ra_lim = sorted([radec0.getRa().asDegrees(), radec1.getRa().asDegrees()])
    dec_lim = sorted([radec0.getDec().asDegrees(), radec1.getDec().asDegrees()])
    injection_catalog = generate_injection_catalog(
        ra_lim=ra_lim,
        dec_lim=dec_lim,
        wcs=wcs,
        number=10,
        source_type="DeltaFunction",
        mag=[10.0, 15.0, 20.0],
    )
    return injection_catalog


def make_test_reference_pipeline():
    """Make a test reference pipeline containing initial single-frame tasks."""
    reference_pipeline = Pipeline("reference_pipeline")
    reference_pipeline.addTask(IsrTask, "isr")
    reference_pipeline.addTask(CharacterizeImageTask, "characterizeImage")
    reference_pipeline.addTask(CalibrateTask, "calibrate")
    reference_pipeline._pipelineIR.labeled_subsets["test_subset"] = LabeledSubset("test_subset", set(), None)
    reference_pipeline.addLabelToSubset("test_subset", "isr")
    reference_pipeline.addLabelToSubset("test_subset", "characterizeImage")
    reference_pipeline.addLabelToSubset("test_subset", "calibrate")
    return reference_pipeline
