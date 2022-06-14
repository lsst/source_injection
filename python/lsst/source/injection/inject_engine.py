# This file is part of pipe tasks
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Engine for injecting sources.
"""
import galsim
import numpy as np

import lsst.geom as geom
from lsst.pex.exceptions import LogicError, InvalidParameterError
from lsst.geom import SpherePoint, radians, Point2D


def inject_sources(exposure, objects, calibFluxRadius=12.0, logger=None):
    """Inject sources into given exposure

    Parameters
    ----------
    exposure : `lsst.afw.image.exposure.exposure.ExposureF`
        The exposure into which the sources should be injected.
    objects : `typing.Iterator` [`tuple` [`lsst.geom.SpherePoint`,
                                          `galsim.GSObject`]]
        An iterator of tuples that contains (or generates) locations and object
        surface brightness profiles to inject.
    calibFluxRadius : `float`, optional
        Aperture radius (in pixels) used to define the calibration for this
        exposure+catalog. This is used to produce the correct instrumental
        fluxes within the radius. The value should match that of the field
        defined in slot_CalibFlux_instFlux.
    logger : `lsst.log.log.log.Log` or `logging.Logger`, optional
        Logger.
    """
    exposure.mask.addMaskPlane("SI")
    bitmask = exposure.mask.getPlaneBitMask("SI")
    if logger:
        logger.info(f"Adding mask plane with bitmask {bitmask}")

    wcs = exposure.getWcs()
    psf = exposure.getPsf()

    bbox = exposure.getBBox()
    fullBounds = galsim.BoundsI(bbox.minX, bbox.maxX, bbox.minY, bbox.maxY)
    gsImg = galsim.Image(exposure.image.array, bounds=fullBounds)

    pixScale = wcs.getPixelScale(bbox.getCenter()).asArcseconds()

    for spt, gsObj in objects:
        pt = wcs.skyToPixel(spt)
        posd = galsim.PositionD(pt.x, pt.y)
        posi = galsim.PositionI(pt.x//1, pt.y//1)
        if logger:
            logger.debug(f"Adding synthetic source at {pt}")

        mat = wcs.linearizePixelToSky(spt, geom.arcseconds).getMatrix()
        gsWCS = galsim.JacobianWCS(mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1])

        # This check is here because sometimes the WCS
        # is multivalued and objects that should not be
        # were being included.
        gsPixScale = np.sqrt(gsWCS.pixelArea())
        if gsPixScale < pixScale/2 or gsPixScale > pixScale*2:
            continue

        try:
            psfArr = psf.computeKernelImage(pt).array
        except InvalidParameterError:
            # Try mapping to nearest point contained in bbox.
            contained_pt = Point2D(
                np.clip(pt.x, bbox.minX, bbox.maxX),
                np.clip(pt.y, bbox.minY, bbox.maxY)
            )
            if pt == contained_pt:  # no difference, so skip immediately
                if logger:
                    logger.info(f"Cannot compute PSF for object at {pt}; skipping")
                continue
            # otherwise, try again with new point
            try:
                psfArr = psf.computeKernelImage(contained_pt).array
            except InvalidParameterError:
                if logger:
                    logger.info(f"Cannot compute PSF for object at {pt}; skipping")
                continue

        apCorr = psf.computeApertureFlux(calibFluxRadius, psf.getAveragePosition())
        psfArr /= apCorr
        gsPSF = galsim.InterpolatedImage(galsim.Image(psfArr), wcs=gsWCS)

        conv = galsim.Convolve(gsObj, gsPSF)
        stampSize = conv.getGoodImageSize(gsWCS.minLinearScale())
        subBounds = galsim.BoundsI(posi).withBorder(stampSize//2)
        subBounds &= fullBounds

        if subBounds.area() > 0:
            subImg = gsImg[subBounds]
            offset = posd - subBounds.true_center
            # Note, for calexp injection, pixel is already part of the PSF and
            # for coadd injection, it's incorrect to include the output pixel.
            # So for both cases, we draw using method='no_pixel'.

            conv.drawImage(
                subImg,
                add_to_image=True,
                offset=offset,
                wcs=gsWCS,
                method='no_pixel'
            )

            subBox = geom.Box2I(
                geom.Point2I(subBounds.xmin, subBounds.ymin),
                geom.Point2I(subBounds.xmax, subBounds.ymax)
            )
            exposure[subBox].mask.array |= bitmask


def _generateGSObjectsFromCatalog(exposure, siCat, galCheckVal, starCheckVal, sourceTypeColumn, logger):
    """Process catalog to generate `galsim.GSObject`s.

    Parameters
    ----------
    exposure : `lsst.afw.image.exposure.exposure.ExposureF`
        The exposure into which the synthetic sources should be added
    siCat : `pandas.core.frame.DataFrame`
        The catalog of synthetic sources to be injected
    galCheckVal : `str`, `bytes` or `int`
        The value set in the sourceType column to specifiy a galaxy object.
    starCheckVal : `str`, `bytes` or `int`
        The value set in the sourceType column to specifiy a star object.

    Yields
    ------
    gsObjects : `generator`
        A generator of tuples of `lsst.geom.SpherePoint` and `galsim.GSObject`.
    """
    wcs = exposure.getWcs()
    photoCalib = exposure.getPhotoCalib()

    logger.info(f"Generating {len(siCat)} sources for injection")

    for (index, row) in siCat.iterrows():
        ra = row["ra"]
        dec = row["dec"]
        skyCoord = SpherePoint(ra, dec, radians)
        xy = wcs.skyToPixel(skyCoord)

        try:
            flux = photoCalib.magnitudeToInstFlux(row["mag"], xy)
        except LogicError:
            continue

        sourceType = row[sourceTypeColumn]
        if sourceType == galCheckVal:
            # GalSim convention: HLR = sqrt(a * b) = a * sqrt(b / a)
            bulge_gs_HLR = row["bulge_semimajor"]*np.sqrt(row["bulge_axis_ratio"])
            bulge = galsim.Sersic(n=row["bulge_n"], half_light_radius=bulge_gs_HLR)
            bulge = bulge.shear(q=row["bulge_axis_ratio"], beta=((90 - row["bulge_pa"])*galsim.degrees))

            disk_gs_HLR = row["disk_semimajor"]*np.sqrt(row["disk_axis_ratio"])
            disk = galsim.Sersic(n=row["disk_n"], half_light_radius=disk_gs_HLR)
            disk = disk.shear(q=row["disk_axis_ratio"], beta=((90 - row["disk_pa"])*galsim.degrees))

            gal = bulge*row["bulge_disk_flux_ratio"] + disk
            gal = gal.withFlux(flux)

            yield skyCoord, gal
        elif sourceType == starCheckVal:
            star = galsim.DeltaFunction()
            star = star.withFlux(flux)
            yield skyCoord, star
        else:
            raise TypeError(f"Unknown sourceType {sourceType}")


# def _generateGSObjectsFromImages(exposure, siCat, logger):
#     """Process catalog to generate `galsim.GSObject` s.

#     Parameters
#     ----------
#     exposure : `lsst.afw.image.exposure.exposure.ExposureF`
#         The exposure into which the synthetic sources should be injected
#     siCat : `pandas.core.frame.DataFrame`
#         The catalog of synthetic sources to be injected

#     Yields
#     ------
#     gsObjects : `generator`
#         A generator of tuples of `lsst.geom.SpherePoint`, `galsim.GSObject`.
#     """
#     band = exposure.getFilterLabel().bandLabel
#     wcs = exposure.getWcs()
#     photoCalib = exposure.getPhotoCalib()

#     logger.info("Processing %d synthetic images", len(siCat))

#     for (index, row) in siCat.iterrows():
#         ra = row["ra"]
#         dec = row["dec"]
#         skyCoord = SpherePoint(ra, dec, radians)
#         xy = wcs.skyToPixel(skyCoord)

#         try:
#             flux = photoCalib.magnitudeToInstFlux(row["mag"], xy)
#         except LogicError:
#             continue

#         imFile = row[band+"imFilename"]
#         try:
#             imFile = imFile.decode("utf-8")
#         except AttributeError:
#             pass
#         imFile = imFile.strip()
#         im = galsim.fits.read(imFile, read_header=True)

#         if self.config.fits_alignment == "wcs":
#             # galsim.fits.read will always attach a WCS to its output. If it
#             # can't find a WCS in the FITS header, then it defaults to
#             # scale = 1.0 arcsec / pix.  So if that's the scale, then we
#             # need to check if it was explicitly set or if it's just the
#             # default.  If it's just the default then we should raise an
#             # exception.
#             if _isWCSGalsimDefault(im.wcs, im.header):
#                 raise RuntimeError(
#                     f"Cannot find WCS in input FITS file {imFile}"
#                 )
#         elif self.config.fits_alignment == "pixel":
#             # Here we need to set im.wcs to the local WCS at the target
#             # position.
#             linWcs = wcs.linearizePixelToSky(skyCoord, geom.arcseconds)
#             mat = linWcs.getMatrix()
#             im.wcs = galsim.JacobianWCS(
#                 mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1]
#             )
#         else:
#             raise ValueError(
#                 f"Unknown fits_alignment type {self.config.fits_alignment}"
#             )

#         obj = galsim.InterpolatedImage(im, calculate_stepk=False)
#         obj = obj.withFlux(flux)
#         yield skyCoord, obj


def addPixCoords(siCat, image):

    """Add pixel coordinates to the catalog of synthetic sources, using RA
    and Dec as an input.

    Parameters
    ----------
    siCat : `pandas.core.frame.DataFrame`
                The catalog of synthetic sources to be injected
    image : `lsst.afw.image.exposure.exposure.ExposureF`
                The image into which the synthetic sources should be injected

    Returns
    -------
    siCat : `pandas.core.frame.DataFrame`
    """
    wcs = image.getWcs()
    ras = siCat["ra"].values
    decs = siCat["dec"].values
    xs, ys = wcs.skyToPixelArray(ras, decs)
    siCat["x"] = xs
    siCat["y"] = ys

    return siCat
