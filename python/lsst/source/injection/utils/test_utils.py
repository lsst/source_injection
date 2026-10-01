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

import astropy.units as u
import numpy as np

from lsst.afw.geom import makeCdMatrix, makeSkyWcs
from lsst.afw.image import ExposureF, makePhotoCalibFromCalibZeroPoint
from lsst.geom import Box2I, Extent2I, Point2D, Point2I, SpherePoint, degrees
from lsst.images import YX, Box, Image, Mask, SkyProjection, TractFrame
from lsst.images.cells import (
    CellCoadd,
    CellGrid,
    CellGridBounds,
    CellIJ,
    CellPointSpreadFunction,
    CoaddProvenance,
    PatchDefinition,
)
from lsst.ip.isr.isrTask import IsrTask
from lsst.meas.algorithms.testUtils import plantSources
from lsst.pipe.base import Pipeline
from lsst.pipe.base.pipelineIR import ContractIR, LabeledSubset
from lsst.pipe.tasks.calibrate import CalibrateTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
from lsst.source.injection import generate_injection_catalog


def make_test_cell_coadd(
    exposure: ExposureF,
    cell_shape: YX,
    psf_shape: tuple[int, int],
    missing: set[CellIJ] | None = None,
    image_unit: u.Unit = u.nJy,
    **kwargs,
) -> CellCoadd:
    """Make a cell coadd out of an exposure.

    Parameters
    ----------
    exposure : `lsst.afw.image.ExposureF`
        An exposure with a simple PSF.
    cell_shape : `lsst.images.YX`
        The shape of each cell.
    psf_shape : `tuple[int, int]`
        The shape of the PSF image array.
    missing : `set[CellIJ]` or `None`
        The set of cells that have invalid PSFs, which will be set to nan.
    image_unit : `astropy.units.Unit`
        The unit of the exposure's image plane.
    kwargs
        Additional keyword arguments to pass to the CellCoadd constructor.

    Returns
    -------
    cell_coadd : `lsst.images.cells.CellCoadd`
        The exposure as a CellCoadd.
    """
    if missing is None:
        missing = set()

    exposure_psf = exposure.psf
    box = Box.from_legacy(exposure.getBBox())
    cell_grid = CellGrid(bbox=box, cell_shape=cell_shape)
    cell_grid_bounds = CellGridBounds(grid=cell_grid, bbox=box, missing=missing)

    # Make the cell PSF grid by evaluating the PSF at the center of each cell
    cell_psfs = []
    for cell_i in range(cell_grid.grid_size.i):
        row_psfs = []
        for cell_j in range(cell_grid.grid_size.j):
            if (cell_ij := CellIJ(cell_i, cell_j)) in missing:
                psf = np.full(psf_shape, np.nan)
            else:
                center = cell_grid.bbox_of(cell_ij).to_legacy().getCenter()
                psf = exposure_psf.computeKernelImage(center).array
                if psf.shape != psf_shape:
                    raise ValueError(f"Exposure PSF image at {center=} has shape={psf.shape} != {psf_shape=}")
            row_psfs.append(psf)
        cell_psfs.append(row_psfs)

    cell_psfs = CellPointSpreadFunction(bounds=cell_grid_bounds, array=np.array(cell_psfs))
    # CellCoadds only take TractFrame for now
    sky_projection = SkyProjection.from_legacy(
        exposure.wcs,
        TractFrame(skymap="dummy", tract=0, bbox=box),
        pixel_bounds=box,
    )
    coadd_contributions = CoaddProvenance.make_empty_contribution_table(n_rows=1)
    coadd_inputs = CoaddProvenance.make_empty_contribution_table(n_rows=1)
    patch = PatchDefinition(
        id=0,
        index=YX(0, 0),
        inner_bbox=box,
        cells=cell_grid,
    )

    cell_coadd = CellCoadd(
        image=Image.from_legacy(exposure.image, unit=image_unit),
        variance=Image.from_legacy(exposure.variance),
        mask=Mask.from_legacy(exposure.mask),
        patch=patch,
        provenance=CoaddProvenance(inputs=coadd_inputs, contributions=coadd_contributions),
        psf=cell_psfs,
        sky_projection=sky_projection,
        **kwargs,
    )
    return cell_coadd


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
    isr_out = "isr.connections.ConnectionsClass(config=isr).outputExposure.name"
    char_inp = "characterizeImage.connections.ConnectionsClass(config=characterizeImage).exposure.name"
    char_out = "characterizeImage.connections.ConnectionsClass(config=characterizeImage).characterized.name"
    calib_inp = "calibrate.connections.ConnectionsClass(config=calibrate).exposure.name"
    # When injecting into a post_isr_image, contract1 should be violated.
    # The make_injection_pipeline utility should be robust to this.
    contract1 = ContractIR(
        contract=f"{isr_out} == {char_inp}",
        msg="isr output == characterizeImage input",
    )
    contract2 = ContractIR(
        contract=f"{char_out} == {calib_inp}",
        msg="characterizeImage output == calibrate input",
    )
    reference_pipeline._pipelineIR.contracts = [contract1, contract2]
    return reference_pipeline
