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

__all__ = ["generate_galsim_objects", "inject_galsim_objects_into_exposure"]

import os
from collections import Counter
from collections.abc import Generator
from typing import Any

import galsim
import numpy as np
import numpy.ma as ma
from astropy.io import fits
from astropy.table import Table
from galsim import GalSimFFTSizeError

from lsst.afw.geom import SkyWcs
from lsst.afw.image import ExposureF, PhotoCalib
from lsst.geom import Box2I, Point2D, Point2I, SpherePoint, arcseconds, degrees
from lsst.pex.exceptions import InvalidParameterError, LogicError


def get_object_data(source_data: dict[str, Any], object_class: galsim.GSObject) -> dict[str, Any]:
    """Assemble a dictionary of allowed keyword arguments and their
    corresponding values to use when constructing a GSObject.

    Parameters
    ----------
    source_data : `dict` [`str`, `Any`]
        Dictionary of source data.
    object_class : `galsim.gsobject.GSObject`
        Class of GSObject to match against.

    Returns
    -------
    object_data : `dict` [`str`, `Any`]
        Dictionary of source data to pass to the GSObject constructor.
    """
    req = getattr(object_class, "_req_params", {})
    opt = getattr(object_class, "_opt_params", {})
    single = getattr(object_class, "_single_params", {})

    object_data = {}

    # Check required args.
    for key in req:
        if key not in source_data:
            raise ValueError(f"Required parameter {key} not found in input catalog.")
        object_data[key] = source_data[key]

    # Optional args.
    for key in opt:
        if key in source_data:
            object_data[key] = source_data[key]

    # Single args.
    for s in single:
        count = 0
        for key in s:
            if key in source_data:
                count += 1
                if count > 1:
                    raise ValueError(f"Only one of {s.keys()} allowed for type {object_class}.")
                object_data[key] = source_data[key]
        if count == 0:
            raise ValueError(f"One of the args {s.keys()} is required for type {object_class}.")

    return object_data


def get_shear_data(
    source_data: dict[str, Any],
    shear_attributes: list[str] = [
        "g1",
        "g2",
        "g",
        "e1",
        "e2",
        "e",
        "eta1",
        "eta2",
        "eta",
        "q",
        "beta",
        "shear",
    ],
) -> dict[str, Any]:
    """Assemble a dictionary of allowed keyword arguments and their
    corresponding values to use when constructing a Shear.

    Parameters
    ----------
    source_data : `dict` [`str`, `Any`]
        Dictionary of source data.
    shear_attributes : `list` [`str`], optional
        List of allowed shear attributes.

    Returns
    -------
    shear_data : `dict` [`str`, `Any`]
        Dictionary of source data to pass to the Shear constructor.
    """
    shear_params = set(shear_attributes) & set(source_data.keys())
    shear_data = {}
    for shear_param in shear_params:
        if shear_param == "beta":
            shear_data.update({shear_param: source_data[shear_param] * galsim.degrees})
        else:
            shear_data.update({shear_param: source_data[shear_param]})
    return shear_data


def generate_galsim_objects(
    injection_catalog: Table,
    photo_calib: PhotoCalib,
    wcs: SkyWcs,
    fits_alignment: str,
    stamp_prefix: str = "",
    logger: Any | None = None,
) -> Generator[tuple[SpherePoint, Point2D, int, galsim.gsobject.GSObject], None, None]:
    """Generate GalSim objects from an injection catalog.

    Parameters
    ----------
    injection_catalog : `astropy.table.Table`
        Table of sources to be injected.
    photo_calib : `lsst.afw.image.PhotoCalib`
        Photometric calibration used to calibrate injected sources.
    wcs : `lsst.afw.geom.SkyWcs`
        WCS used to calibrate injected sources.
    fits_alignment : `str`
        Alignment of the FITS image to the WCS. Allowed values: "wcs", "pixel".
    stamp_prefix : `str`
        Prefix to add to the stamp name.
    logger : `lsst.utils.logging.LsstLogAdapter`, optional
        Logger to use for logging messages.

    Yields
    ------
    sky_coords : `lsst.geom.SpherePoint`
        RA/Dec coordinates of the source.
    pixel_coords : `lsst.geom.Point2D`
        Pixel coordinates of the source.
    draw_size : `int`
        Size of the stamp to draw.
    object : `galsim.gsobject.GSObject`
        A fully specified and transformed GalSim object.
    """
    if logger:
        source_types = Counter(injection_catalog["source_type"])  # type: ignore
        grammar0 = "source" if len(injection_catalog) == 1 else "sources"
        grammar1 = "type" if len(source_types) == 1 else "types"
        logger.info(
            "Generating %d injection %s consisting of %d unique %s: %s.",
            len(injection_catalog),
            grammar0,
            len(source_types),
            grammar1,
            ", ".join(f"{k}({v})" for k, v in source_types.items()),
        )
    for source_data_full in injection_catalog:
        items = dict(source_data_full).items()  # type: ignore
        source_data = {k: v for (k, v) in items if v is not ma.masked}
        try:
            sky_coords = SpherePoint(float(source_data["ra"]), float(source_data["dec"]), degrees)
        except KeyError:
            sky_coords = wcs.pixelToSky(float(source_data["x"]), float(source_data["y"]))
        try:
            pixel_coords = Point2D(source_data["x"], source_data["y"])
        except KeyError:
            pixel_coords = wcs.skyToPixel(sky_coords)
        try:
            inst_flux = photo_calib.magnitudeToInstFlux(source_data["mag"], pixel_coords)
        except LogicError:
            continue
        try:
            draw_size = int(source_data["draw_size"])
        except KeyError:
            draw_size = 0

        if source_data["source_type"] == "Stamp":
            stamp_file = stamp_prefix + source_data["stamp"]
            object = make_galsim_stamp(stamp_file, fits_alignment, wcs, sky_coords, inst_flux)
        elif source_data["source_type"] == "Trail":
            object = make_galsim_trail(source_data, wcs, sky_coords, inst_flux)
        else:
            object = make_galsim_object(source_data, source_data["source_type"], inst_flux, logger)

        yield sky_coords, pixel_coords, draw_size, object


def make_galsim_object(
    source_data: dict[str, Any],
    source_type: str,
    inst_flux: float,
    logger: Any | None = None,
) -> galsim.gsobject.GSObject:
    """Make a generic GalSim object from a collection of source data.

    Parameters
    ----------
    source_data : `dict` [`str`, `Any`]
        Dictionary of source data.
    source_type : `str`
        Type of the source, corresponding to a GalSim class.
    inst_flux : `float`
        Instrumental flux of the source.
    logger : `~lsst.utils.logging.LsstLogAdapter`, optional

    Returns
    -------
    object : `galsim.gsobject.GSObject`
        A fully specified and transformed GalSim object.
    """
    # Populate the non-shear and non-flux parameters.
    object_class = getattr(galsim, source_type)
    object_data = get_object_data(source_data, object_class)
    object = object_class(**object_data)
    # Create a version of the object with an area-preserving shear applied.
    shear_data = get_shear_data(source_data)
    if shear_data:
        try:
            object = object.shear(**shear_data)
        except TypeError as err:
            if logger:
                logger.warning("Cannot apply shear to GalSim object: %s", err)
            pass
    # Apply the instrumental flux and return.
    object = object.withFlux(inst_flux)
    return object


def make_galsim_trail(
    source_data: dict[str, Any],
    wcs: SkyWcs,
    sky_coords: SpherePoint,
    inst_flux: float,
    trail_thickness: float = 1e-6,
) -> galsim.gsobject.GSObject:
    """Make a trail with GalSim from a collection of source data.

    Parameters
    ----------
    source_data : `dict` [`str`, `Any`]
        Dictionary of source data.
    wcs : `lsst.afw.geom.SkyWcs`
        World coordinate system.
    sky_coords : `lsst.geom.SpherePoint`
        Sky coordinates of the source.
    inst_flux : `float`
        Instrumental flux of the source.
    trail_thickness : `float`
        Thickness of the trail in pixels.

    Returns
    -------
    object : `galsim.gsobject.GSObject`
        A fully specified and transformed GalSim object.
    """
    # Make a 'thin' box to mimic a line surface brightness profile of default
    # thickness = 1e-6 (i.e., much thinner than a pixel)
    object = galsim.Box(source_data["trail_length"], trail_thickness)
    try:
        object = object.rotate(source_data["beta"] * galsim.degrees)
    except KeyError:
        pass
    object = object.withFlux(inst_flux * source_data["trail_length"])  # type: ignore
    # GalSim objects are assumed to be in sky-coords. As we want the trail to
    # appear as defined above in image-coords, we must transform it here.
    linear_wcs = wcs.linearizePixelToSky(sky_coords, arcseconds)
    mat = linear_wcs.getMatrix()
    object = object.transform(mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1])
    # Rescale flux; transformation preserves surface brightness, not flux.
    object *= 1.0 / np.abs(mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0])
    return object  # type: ignore


def make_galsim_stamp(
    stamp_file: str,
    fits_alignment: str,
    wcs: SkyWcs,
    sky_coords: SpherePoint,
    inst_flux: float,
) -> galsim.gsobject.GSObject:
    """Make a postage stamp with GalSim from a FITS file.

    Parameters
    ----------
    stamp_file : `str`
        Path to the FITS file containing the postage stamp.
    fits_alignment : `str`
        Alignment of the FITS image to the WCS. Allowed values: "wcs", "pixel".
    wcs : `lsst.afw.geom.SkyWcs`
        World coordinate system.
    sky_coords : `lsst.geom.SpherePoint`
        Sky coordinates of the source.
    inst_flux : `float`
        Instrumental flux of the source.

    Returns
    -------
    object: `galsim.gsobject.GSObject`
        A fully specified and transformed GalSim object.
    """
    stamp_file = stamp_file.strip()
    if os.path.exists(stamp_file):
        with fits.open(stamp_file) as hdul:
            hdu_images = [hdu.is_image and hdu.size > 0 for hdu in hdul]  # type: ignore
        if any(hdu_images):
            stamp_data = galsim.fits.read(stamp_file, read_header=True, hdu=np.where(hdu_images)[0][0])
        else:
            raise RuntimeError(f"Cannot find image in input FITS file {stamp_file}.")
        match fits_alignment:
            case "wcs":
                # galsim.fits.read will always attach a WCS to its output.
                # If it can't find a WCS in the FITS header, then it
                # defaults to scale = 1.0 arcsec / pix. If that's the scale
                # then we need to check if it was explicitly set or if it's
                # just the default. If it's just the default then we should
                # raise an exception.
                if is_wcs_galsim_default(stamp_data.wcs, stamp_data.header):  # type: ignore
                    raise RuntimeError(f"Cannot find WCS in input FITS file {stamp_file}")
            case "pixel":
                # We need to set stamp_data.wcs to the local target WCS.
                linear_wcs = wcs.linearizePixelToSky(sky_coords, arcseconds)
                mat = linear_wcs.getMatrix()
                stamp_data.wcs = galsim.JacobianWCS(  # type: ignore
                    mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1]
                )
        object = galsim.InterpolatedImage(stamp_data, calculate_stepk=False)
        object = object.withFlux(inst_flux)
        return object  # type: ignore
    else:
        raise RuntimeError(f"Cannot locate input FITS postage stamp {stamp_file}.")


def is_wcs_galsim_default(
    wcs: galsim.fitswcs.GSFitsWCS,
    hdr: galsim.fits.FitsHeader,
) -> bool:
    """Decide if wcs = galsim.PixelScale(1.0) is explicitly present in header,
    or if it's just the GalSim default.

    Parameters
    ----------
    wcs : galsim.fitswcs.GSFitsWCS
        Potentially default WCS.
    hdr : galsim.fits.FitsHeader
        Header as read in by GalSim.

    Returns
    -------
    is_default : bool
        True if default, False if explicitly set in header.
    """
    if wcs != galsim.PixelScale(1.0):
        return False
    if hdr.get("GS_WCS") is not None:
        return False
    if hdr.get("CTYPE1", "LINEAR") == "LINEAR":
        return not any(k in hdr for k in ["CD1_1", "CDELT1"])
    for wcs_type in galsim.fitswcs.fits_wcs_types:
        # If one of these succeeds, then assume result is explicit.
        try:
            wcs_type._readHeader(hdr)
            return False
        except Exception:
            pass
    else:
        return not any(k in hdr for k in ["CD1_1", "CDELT1"])


def get_gains_from_amps(
    exposure: ExposureF,
) -> tuple[list[Box2I], list[float]] | tuple[None, None]:
    """Retrieve the bounding boxes and gains for each amplifier, if applicable.
    If this cannot be done, return None.

    Parameters
    ----------
    exposure : `lsst.afw.image.ExposureF`
        The exposure to inject synthetic sources into.

    Returns
    -------
    bboxes : `list` [`lsst.geom.Box2I`]
        Bounding boxes for each amplifier region in the exposure.
    gains : `list` [`float`]
        The gain recorded for each amplifier in this exposure.
    """
    try:
        amplifiers = exposure.getDetector().getAmplifiers()
        bboxes = [amplifier.getBBox() for amplifier in amplifiers]
        gains = [amplifier.getGain() for amplifier in amplifiers]
    except AttributeError:
        return None, None
    return bboxes, gains


def infer_gain_from_image(
    exposure: ExposureF,
    bad_mask_names: list[str] | None = None,
) -> float:
    """Fit a straight line to image flux vs. estimated variance in each pixel
    of an exposure, and from that infer an estimated average gain.

    Parameters
    ----------
    exposure : `lsst.afw.image.ExposureF`
        The exposure to inject synthetic sources into.
    bad_mask_names : `list` [`str`], optional
        Names of any mask plane bits that could indicate a pathological flux
        vs. variance relationship.

    Returns
    -------
    gain : `float`
        The estimated average gain across all pixels.
    """
    image_arr = exposure.image.array
    var_arr = exposure.variance.array
    mask_arr = exposure.mask.array
    good = np.isfinite(image_arr) & np.isfinite(var_arr)
    if bad_mask_names is not None:
        bad_mask_bit_mask = exposure.mask.getPlaneBitMask(bad_mask_names)
        good &= (mask_arr.astype(int) & bad_mask_bit_mask) == 0
    fit = np.polyfit(image_arr[good], var_arr[good], deg=1)
    return 1.0 / fit[0]


def get_gain_map(
    exposure: ExposureF,
    full_bounds: galsim.BoundsI,
    bad_mask_names: list[str] | None = None,
) -> tuple[galsim.Image, bool]:
    """Retrieve gains for each pixel.

    Parameters
    ----------
    exposure : `lsst.afw.image.ExposureF`
        The exposure to inject synthetic sources into.
    full_bounds : `galsim.BoundsI`
        Bounding box covering the entire exposure.
    bad_mask_names : `list` [`str`], optional
        Names of any mask plane bits that could indicate a pathological flux
        vs. variance relationship.

    Returns
    -------
    gain_map : `galsim.Image`
        The gain to use in each pixel.
    is_single_CCD : `bool`
        Whether the exposure represents a single CCD, having well-defined
        amplifier regions with recorded gains.
    """
    gain_map = galsim.Image(full_bounds, dtype=float)
    is_single_CCD = True

    bboxes, gains = get_gains_from_amps(exposure)
    if bboxes is None:
        # Separate amplifier bounding boxes not found.
        # This exposure does not represent a single CCD.
        is_single_CCD = False
        bboxes = [full_bounds]
        gains = [infer_gain_from_image(exposure, bad_mask_names)]

    assert gains is not None  # Make mypy happy.
    for bbox, gain in zip(bboxes, gains):
        if is_single_CCD:
            bounds = galsim.BoundsI(bbox.minX, bbox.maxX, bbox.minY, bbox.maxY)
        else:
            bounds = galsim.BoundsI(bbox.getXMin(), bbox.getXMax(), bbox.getYMin(), bbox.getYMax())
        gain_map[bounds].array[:] = 1.0 / gain
    return gain_map, is_single_CCD


def add_noise_to_galsim_image(
    image_template: galsim.Image,
    var_template: galsim.Image,
    is_single_CCD: bool,
    gain_map: galsim.Image,
    object_common_bounds: galsim.BoundsI,
    noise_seed: int | None = None,
    logger: Any | None = None,
) -> tuple[galsim.Image, galsim.Image]:
    """Add noise to a supplied image, and adjust the estimated variance
    accordingly. If the image represents a single CCD exposure, add Poisson
    shot noise. Otherwise add Gaussian noise.

    Parameters
    ----------
    image_template : `galsim.Image`
        The image to add noise to.
    var_template : `galsim.Image`
        The variance to use for the noise generating process in each pixel.
    is_single_CCD : `bool`
        Whether the exposure represents a single CCD, having well-defined
        amplifier regions with recorded gains.
    gain_map : `galsim.Image`
        The gain to use in each pixel of the full exposure.
    object_common_bounds : `galsim.BoundsI`
        Specific region of the full exposure to which this image belongs.
    noise_seed : `int`, optional
        Random number generator seed to use for the noise realization.
    logger : `lsst.utils.logging.LsstLogAdapter`, optional
        Logger to use for logging messages.

    Returns
    -------
    image_template : `galsim.Image`
        The image with noise added.
    var_template : `galsim.Image`
        The estimated variance, adjusted in reponse to the added noise.
    """
    noise_template = var_template.copy()
    negative_variance = var_template.array < 0
    nonfinite_variance = ~np.isfinite(var_template.array)
    if np.any(negative_variance):
        if logger:
            logger.debug("Setting negative-variance pixels to 0 variance for noise generation.")
        noise_template.array[negative_variance] = 0
    if np.any(nonfinite_variance):
        if logger:
            logger.debug("Setting non-finite-variance pixels to 0 variance for noise generation.")
        noise_template.array[nonfinite_variance] = 0

    rng = galsim.BaseDeviate(noise_seed)

    if is_single_CCD:
        # Treat noise_template, which is the expected amount of
        # injected flux divided by the gain, as the Poisson mean of
        # the number of photons collected by each pixel.
        noise = galsim.PoissonNoise(rng)
        noise_template.addNoise(noise)
        # Scale the photon count by the gain to get the amount of
        # image flux to inject.
        image_template = noise_template.copy()
        image_template /= gain_map[object_common_bounds]
        # Restore the original variance values
        # of the bad-variance pixels.
        bad_var = negative_variance | nonfinite_variance
        noise_template.array[bad_var] = var_template.array[bad_var]
        var_template = noise_template
    else:
        noise = galsim.VariableGaussianNoise(rng, noise_template)
        image_template.addNoise(noise)
        var_template = image_template.copy()
        var_template *= gain_map[object_common_bounds]

    return image_template, var_template


def inject_galsim_objects_into_exposure(
    exposure: ExposureF,
    objects: Generator[tuple[SpherePoint, Point2D, int, galsim.gsobject.GSObject], None, None],
    mask_plane_name: str = "INJECTED",
    calib_flux_radius: float = 12.0,
    draw_size_max: int = 1000,
    inject_variance: bool = True,
    add_noise: bool = True,
    noise_seed: int = 0,
    bad_mask_names: list[str] | None = None,
    logger: Any | None = None,
) -> tuple[list[int], list[galsim.BoundsI], list[bool], list[bool]]:
    """Inject sources into given exposure using GalSim.

    Parameters
    ----------
    exposure : `lsst.afw.image.ExposureF`
        The exposure to inject synthetic sources into.
    objects : `Generator` [`tuple`, None, None]
        An iterator of tuples that contains (or generates) locations and object
        surface brightness profiles to inject. The tuples should contain the
        following elements: `lsst.geom.SpherePoint`, `lsst.geom.Point2D`,
        `int`, `galsim.gsobject.GSObject`.
    mask_plane_name : `str`
        Name of the mask plane to use for the injected sources.
    calib_flux_radius : `float`
        Radius in pixels to use for the flux calibration. This is used to
        produce the correct instrumental fluxes within the radius. The value
        should match that of the field defined in slot_CalibFlux_instFlux.
    draw_size_max : `int`
        Maximum allowed size of the drawn object. If the object is larger than
        this, the draw size will be clipped to this size.
    inject_variance : `bool`
        Whether, when injecting flux into the image plane, to inject a
        corresponding amount of variance into the variance plane.
    add_noise : `bool`
        Whether to randomly vary the amount of injected image flux in each
        pixel by an amount consistent with the amount of injected variance.
    noise_seed : `int`
        Seed for generating the noise of the first injected object. The seed
        actually used increases by 1 for each subsequent object, to ensure
        independent noise realizations.
    bad_mask_names : `list[str]`, optional
        List of mask plane names indicating pixels to ignore when fitting flux
        vs variance in preparation for variance plane modification. If None,
        then the all pixels are used.
    logger : `lsst.utils.logging.LsstLogAdapter`, optional
        Logger to use for logging messages.

    Returns
    -------
    draw_sizes : `list` [`int`]
        Draw sizes of the injected sources.
    common_bounds : `list` [`galsim.BoundsI`]
        Common bounds of the drawn objects.
    fft_size_errors : `list` [`bool`]
        Boolean flags indicating whether a GalSimFFTSizeError was raised.
    psf_compute_errors : `list` [`bool`]
        Boolean flags indicating whether a PSF computation error was raised.
    """
    exposure.mask.addMaskPlane(mask_plane_name)
    mask_plane_core_name = mask_plane_name + "_CORE"
    exposure.mask.addMaskPlane(mask_plane_core_name)
    if logger:
        logger.info(
            "Adding %s and %s mask planes to the exposure.",
            mask_plane_name,
            mask_plane_core_name,
        )
    psf = exposure.getPsf()
    wcs = exposure.getWcs()
    bbox = exposure.getBBox()
    full_bounds = galsim.BoundsI(bbox.minX, bbox.maxX, bbox.minY, bbox.maxY)
    galsim_image = galsim.Image(exposure.image.array, bounds=full_bounds)
    galsim_variance = galsim.Image(exposure.variance.array, bounds=full_bounds)
    pixel_scale = wcs.getPixelScale(bbox.getCenter()).asArcseconds()

    gain_map, is_single_CCD = get_gain_map(exposure, full_bounds, bad_mask_names)

    draw_sizes: list[int] = []
    common_bounds: list[galsim.BoundsI] = []
    fft_size_errors: list[bool] = []
    psf_compute_errors: list[bool] = []
    for i, (sky_coords, pixel_coords, draw_size, object) in enumerate(objects):
        # Instantiate default returns in case of early exit from this loop.
        draw_sizes.append(0)
        common_bounds.append(galsim.BoundsI())
        fft_size_errors.append(False)
        psf_compute_errors.append(False)

        # Get spatial coordinates and WCS.
        posd = galsim.PositionD(pixel_coords.x, pixel_coords.y)
        posi = galsim.PositionI(pixel_coords.x // 1, pixel_coords.y // 1)
        if logger:
            logger.debug(f"Injecting synthetic source at {pixel_coords}.")
        mat = wcs.linearizePixelToSky(sky_coords, arcseconds).getMatrix()
        galsim_wcs = galsim.JacobianWCS(mat[0, 0], mat[0, 1], mat[1, 0], mat[1, 1])

        # This check is here because sometimes the WCS is multivalued and
        # objects that should not be included were being included.
        galsim_pixel_scale = np.sqrt(galsim_wcs.pixelArea())
        if galsim_pixel_scale < pixel_scale / 2 or galsim_pixel_scale > pixel_scale * 2:
            continue

        # Compute the PSF at the object location.
        try:
            psf_array = psf.computeKernelImage(pixel_coords).array
        except InvalidParameterError:
            # Try mapping to nearest point contained in bbox.
            contained_point = Point2D(
                np.clip(pixel_coords.x, bbox.minX, bbox.maxX), np.clip(pixel_coords.y, bbox.minY, bbox.maxY)
            )
            if pixel_coords == contained_point:  # no difference, so skip immediately
                psf_compute_errors[i] = True
                if logger:
                    logger.debug("Cannot compute PSF for object at %s; flagging and skipping.", sky_coords)
                continue
            # Otherwise, try again with new point.
            try:
                psf_array = psf.computeKernelImage(contained_point).array
            except InvalidParameterError:
                psf_compute_errors[i] = True
                if logger:
                    logger.debug("Cannot compute PSF for object at %s; flagging and skipping.", sky_coords)
                continue

        # Compute the aperture corrected PSF interpolated image.
        aperture_correction = psf.computeApertureFlux(calib_flux_radius, psf.getAveragePosition())
        psf_array /= aperture_correction
        galsim_psf = galsim.InterpolatedImage(galsim.Image(psf_array), wcs=galsim_wcs)

        # Convolve the object with the PSF and generate draw size.
        conv = galsim.Convolve(object, galsim_psf)
        if draw_size == 0:
            draw_size = conv.getGoodImageSize(galsim_wcs.minLinearScale())  # type: ignore
        injection_draw_size = int(draw_size)
        injection_core_size = 3
        if draw_size_max > 0 and injection_draw_size > draw_size_max:
            if logger:
                logger.warning(
                    "Clipping draw size for object at %s from %d to %d pixels.",
                    sky_coords,
                    injection_draw_size,
                    draw_size_max,
                )
            injection_draw_size = draw_size_max
        draw_sizes[i] = injection_draw_size
        if injection_core_size > injection_draw_size:
            if logger:
                logger.debug(
                    "Clipping core size for object at %s from %d to %d pixels.",
                    sky_coords,
                    injection_core_size,
                    injection_draw_size,
                )
            injection_core_size = injection_draw_size
        sub_bounds = galsim.BoundsI(posi).withBorder(injection_draw_size // 2)
        object_common_bounds = full_bounds & sub_bounds
        common_bounds[i] = object_common_bounds  # type: ignore

        # Inject the source if there is any overlap.
        if object_common_bounds.area() > 0:
            common_image = galsim_image[object_common_bounds]
            common_variance = galsim_variance[object_common_bounds]

            offset = posd - object_common_bounds.true_center
            # Note, for calexp injection, pixel is already part of the PSF and
            # for coadd injection, it's incorrect to include the output pixel.
            # So for both cases, we draw using method='no_pixel'.
            image_template = common_image.copy()
            try:
                image_template = conv.drawImage(
                    image_template, add_to_image=False, offset=offset, wcs=galsim_wcs, method="no_pixel"
                )
            except GalSimFFTSizeError as err:
                fft_size_errors[i] = True
                if logger:
                    logger.debug(
                        "GalSimFFTSizeError raised for object at %s; flagging and skipping.\n%s",
                        sky_coords,
                        err,
                    )
                continue

            # The amount of additional variance to inject,
            # if inject_variance is true.
            var_template = image_template.copy()
            var_template *= gain_map[object_common_bounds]

            if add_noise:
                image_template, var_template = add_noise_to_galsim_image(
                    image_template,
                    var_template,
                    is_single_CCD,
                    gain_map,
                    object_common_bounds,
                    noise_seed,
                    logger,
                )

            common_image += image_template
            if inject_variance:
                common_variance += var_template

            # Increment the seed so different noise is generated for different
            # objects.
            noise_seed += 1

            common_box = Box2I(
                Point2I(object_common_bounds.xmin, object_common_bounds.ymin),
                Point2I(object_common_bounds.xmax, object_common_bounds.ymax),
            )
            bitvalue = exposure.mask.getPlaneBitMask(mask_plane_name)
            exposure[common_box].mask.array |= bitvalue
            # Add a 3 x 3 pixel mask centered on the object. The mask must be
            # large enough to always identify the core/peak of the injected
            # source, but small enough that it rarely overlaps real sources.
            sub_bounds_core = galsim.BoundsI(posi).withBorder(injection_core_size // 2)
            object_common_bounds_core = full_bounds & sub_bounds_core
            if object_common_bounds_core.area() > 0:
                common_box_core = Box2I(
                    Point2I(object_common_bounds_core.xmin, object_common_bounds_core.ymin),
                    Point2I(object_common_bounds_core.xmax, object_common_bounds_core.ymax),
                )
                bitvalue_core = exposure.mask.getPlaneBitMask(mask_plane_core_name)
                exposure[common_box_core].mask.array |= bitvalue_core
        else:
            if logger:
                logger.debug("No area overlap for object at %s; flagging and skipping.", sky_coords)

    return draw_sizes, common_bounds, fft_size_errors, psf_compute_errors
