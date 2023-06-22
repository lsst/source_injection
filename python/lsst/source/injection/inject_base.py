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

__all__ = ["BaseInjectConnections", "BaseInjectConfig", "BaseInjectTask"]

from typing import cast

import galsim
import numpy as np
import numpy.ma as ma
from astropy import units
from astropy.table import Table, hstack, vstack
from astropy.units import Quantity, UnitConversionError
from lsst.afw.image.exposure.exposureUtils import bbox_contains_sky_coords
from lsst.pex.config import ChoiceField, Field
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base.connectionTypes import PrerequisiteInput

from .inject_engine import generate_galsim_objects, inject_galsim_objects_into_exposure


class BaseInjectConnections(
    PipelineTaskConnections,
    dimensions=("instrument",),
    defaultTemplates={
        "injection_prefix": "injection_",
        "injected_prefix": "injected_",
    },
):
    """Base connections for source injection tasks."""

    injection_catalogs = PrerequisiteInput(
        doc="Set of catalogs of sources to draw inputs from.",
        name="{injection_prefix}catalog",
        dimensions=("htm7", "band"),
        storageClass="ArrowAstropy",
        minimum=0,
        multiple=True,
    )


class BaseInjectConfig(PipelineTaskConfig, pipelineConnections=BaseInjectConnections):
    """Base configuration for source injection tasks."""

    # Catalog manipulation options.
    process_all_data_ids = Field[bool](
        doc="If True, all input data IDs will be processed, even those where no synthetic sources were "
        "identified for injection. In such an eventuality this returns a clone of the input image, renamed "
        "to the *output_exposure* connection name and with an empty *mask_plane_name* mask plane attached.",
        default=False,
    )
    trim_padding = Field[int](
        doc="Size of the pixel padding surrounding the image. Only those synthetic sources with a centroid "
        "falling within the ``image + trim_padding`` region will be considered for source injection.",
        default=100,
        optional=True,
    )
    selection = Field[str](
        doc="A string that can be evaluated as a boolean expression to select rows in the input injection "
        "catalog. To make use of this configuration option, the internal object name ``injection_catalog`` "
        "must be used. For example, to select all sources with a magnitude in the range 20.0 < mag < 25.0, "
        "set ``selection=\"(injection_catalog['mag'] > 20.0) & (injection_catalog['mag'] < 25.0)\"``. "
        "The ``{visit}`` field will be substituted for the current visit ID of the exposure being processed. "
        "For example, to select only visits that match a user-supplied visit column in the input injection "
        "catalog, set ``selection=\"np.isin(injection_catalog['visit'], {visit})\"``.",
        optional=True,
    )

    # General configuration options.
    mask_plane_name = Field[str](
        doc="Name assigned to the injected mask plane which is attached to the output exposure.",
        default="INJECTED",
    )
    calib_flux_radius = Field[float](
        doc="Aperture radius (in pixels) that was used to define the calibration for this image+catalog. "
        "This will be used to produce the correct instrumental fluxes within the radius. "
        "This value should match that of the field defined in ``slot_CalibFlux_instFlux``.",
        default=12.0,
    )
    fits_alignment = ChoiceField[str](  # type: ignore
        doc="How should injections from FITS files be aligned?",
        dtype=str,
        allowed={
            "wcs": (
                "Input image will be transformed such that the local WCS in the FITS header matches the "
                "local WCS in the target image. I.e., North, East, and angular distances in the input image "
                "will match North, East, and angular distances in the target image."
            ),
            "pixel": (
                "Input image will **not** be transformed. Up, right, and pixel distances in the input image "
                "will match up, right and pixel distances in the target image."
            ),
        },
        default="pixel",
    )
    stamp_prefix = Field[str](
        doc="String to prefix to the entries in the *col_stamp* column, for example, a directory path.",
        default="",
    )

    # Custom column names.
    col_ra = Field[str](
        doc="Column name for right ascension (in degrees).",
        default="ra",
    )
    col_dec = Field[str](
        doc="Column name for declination (in degrees).",
        default="dec",
    )
    col_mag = Field[str](
        doc="Column name for magnitude.",
        default="mag",
    )
    col_source_type = Field[str](
        doc="Column name for the source type used in the input catalog. Must match one of the surface "
        "brightness profiles defined by GalSim.",
        default="source_type",
    )
    col_stamp = Field[str](
        doc="Column name to identify FITS file postage stamps for direct injection. The strings in this "
        "column will be prefixed with a string given in *stamp_prefix*, to assist in providing the full "
        "path to a FITS file.",
        default="stamp",
    )
    col_draw_size = Field[str](
        doc="Column name providing pixel size of the region into which the source profile will be drawn. If "
        "this column is not provided as an input, the GalSim method ``getGoodImageSize`` will be used "
        "instead.",
        default="draw_size",
    )
    col_trail_length = Field[str](
        doc="Column name for specifying a satellite trail length (in pixels).",
        default="trail_length",
    )

    def setDefaults(self):
        super().setDefaults()


class BaseInjectTask(PipelineTask):
    """Base class for injecting sources into images."""

    _DefaultName = "baseInjectTask"
    ConfigClass = BaseInjectConfig

    def run(self, injection_catalogs, input_exposure, psf, photo_calib, wcs):
        """Inject sources into an image.

        Parameters
        ----------
        injection_catalogs : `list` [`lsst.daf.butler.DeferredDatasetHandle`]
            Tract level injection catalogs that potentially cover the named
            input exposure.
        input_exposure : `lsst.afw.image.exposure.ExposureF`
            The exposure sources will be injected into.
        psf: `lsst.meas.algorithms.ImagePsf`
            PSF model.
        photo_calib : `lsst.afw.image.photoCalib.PhotoCalib`
            Photometric calibration used to calibrate injected sources.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS used to calibrate injected sources.

        Returns
        -------
        result_struct : `lsst.pipe.base.struct.Struct`
            contains : outputExposure : `lsst.afw.image.exposure.ExposureF`
                       outputCatalog : `lsst.afw.table.source.SourceCatalog`
        """
        self.config = cast(BaseInjectConfig, self.config)

        # Attach potential externally calibrated datasets to input_exposure.
        # Keep originals so we can reset at the end.
        original_psf = input_exposure.getPsf()
        original_photo_calib = input_exposure.getPhotoCalib()
        original_wcs = input_exposure.getWcs()
        input_exposure.setPsf(psf)
        input_exposure.setPhotoCalib(photo_calib)
        input_exposure.setWcs(wcs)

        # Make empty table if none supplied to support process_all_data_ids.
        if len(injection_catalogs) == 0:
            injection_catalogs = [Table()]

        # Consolidate injection catalogs and compose main injection catalog.
        injection_catalog = self._compose_injection_catalog(injection_catalogs)

        # Mapping between standard column names and configured names/units.
        column_mapping = {
            "ra": (self.config.col_ra, units.deg),
            "dec": (self.config.col_dec, units.deg),
            "source_type": (self.config.col_source_type, None),
            "mag": (self.config.col_mag, units.mag),
            "stamp": (self.config.col_stamp, None),
            "draw_size": (self.config.col_draw_size, units.pix),
            "trail_length": (self.config.col_trail_length, units.pix),
        }

        # Standardize injection catalog column names and units.
        injection_catalog = self._standardize_columns(
            injection_catalog,
            column_mapping,
            input_exposure.getWcs().getPixelScale().asArcseconds(),
        )

        # Clean the injection catalog of sources which are not injectable.
        injection_catalog = self._clean_sources(injection_catalog, input_exposure)

        # Injection binary flag lookup dictionary.
        binary_flags = {
            "MAG_BAD": 0,
            "TYPE_UNKNOWN": 1,
            "SERSIC_EXTREME": 2,
            "NO_OVERLAP": 3,
            "FFT_SIZE_ERROR": 4,
            "PSF_COMPUTE_ERROR": 5,
        }

        # Flag sources in the injection catalog prior to source injection.
        injection_catalog = self._flag_sources(injection_catalog, input_exposure, binary_flags)

        # Inject sources into input_exposure.
        good_injections: list[bool] = injection_catalog["injection_flag"] == 0
        good_injections_index = [i for i, val in enumerate(good_injections) if val]
        num_injected_sources = np.sum(good_injections)
        if num_injected_sources > 0:
            object_generator = generate_galsim_objects(
                injection_catalog=injection_catalog[good_injections],
                photo_calib=photo_calib,
                wcs=wcs,
                fits_alignment=self.config.fits_alignment,
                stamp_prefix=self.config.stamp_prefix,
                logger=self.log,
            )
            (
                draw_sizes,
                common_bounds,
                fft_size_errors,
                psf_compute_errors,
            ) = inject_galsim_objects_into_exposure(
                input_exposure,
                object_generator,
                mask_plane_name=self.config.mask_plane_name,
                calib_flux_radius=self.config.calib_flux_radius,
                draw_size_max=10000,  # TODO: replace draw_size logic with GS logic.
                logger=self.log,
            )
            # Add inject_galsim_objects_into_exposure outputs into output cat.
            common_areas = [x.area() if x is not None else None for x in common_bounds]
            for i, (draw_size, common_area, fft_size_error, psf_compute_error) in enumerate(
                zip(draw_sizes, common_areas, fft_size_errors, psf_compute_errors)
            ):
                injection_catalog["injection_draw_size"][good_injections_index[i]] = draw_size
                if common_area == 0:
                    injection_catalog["injection_flag"][good_injections_index[i]] += (
                        2 ** binary_flags["NO_OVERLAP"]
                    )
                if fft_size_error:
                    injection_catalog["injection_flag"][good_injections_index[i]] += (
                        2 ** binary_flags["FFT_SIZE_ERROR"]
                    )
                if psf_compute_error:
                    injection_catalog["injection_flag"][good_injections_index[i]] += (
                        2 ** binary_flags["PSF_COMPUTE_ERROR"]
                    )
        elif num_injected_sources == 0 and self.config.process_all_data_ids:
            self.log.warning("No sources to be injected for this DatasetRef; processing anyway.")
            input_exposure.mask.addMaskPlane(self.config.mask_plane_name)
        else:
            raise RuntimeError(
                "No sources to be injected for this DatasetRef, and process_all_data_ids is False."
            )

        # Restore original input_exposure calibrated data.
        input_exposure.setPsf(original_psf)
        input_exposure.setPhotoCalib(original_photo_calib)
        input_exposure.setWcs(original_wcs)

        # Add injection provenance and injection flags metadata.
        metadata = input_exposure.getMetadata()
        input_dataset_type = self.config.connections.input_exposure.format(**self.config.connections.toDict())
        metadata.set("INJECTED", input_dataset_type, "Initial source injection dataset type")
        for flag, value in sorted(binary_flags.items(), key=lambda item: item[1]):
            injection_catalog.meta[flag] = value

        output_struct = Struct(output_exposure=input_exposure, output_catalog=injection_catalog)
        return output_struct

    def _compose_injection_catalog(self, injection_catalogs):
        """Consolidate injection catalogs and compose main injection catalog.

        If multiple injection catalogs are input, all catalogs are
        concatenated together.

        A running injection_id, specific to this dataset ref, is assigned to
        each source in the output injection catalog if not provided.

        Parameters
        ----------
        injection_catalogs : `list` [`lsst.daf.butler.DeferredDatasetHandle`]
            Set of synthetic source catalogs to concatenate.

        Returns
        -------
        injection_catalog : `astropy.table.Table`
            Catalog of sources to be injected.
        """
        self.config = cast(BaseInjectConfig, self.config)

        # Generate injection IDs (if not provided) and injection flag column.
        injection_data = vstack(injection_catalogs)
        if "injection_id" in injection_data.columns:
            injection_id = injection_data["injection_id"]
            injection_data.remove_column("injection_id")
        else:
            injection_id = range(len(injection_data))
        injection_header = Table(
            dict(
                injection_id=injection_id,
                injection_flag=np.zeros(len(injection_data), dtype=int),
                injection_draw_size=np.zeros(len(injection_data), dtype=int),
            )
        )

        # Construct final injection catalog.
        injection_catalog = hstack([injection_header, injection_data])
        if "source_type" in injection_catalog.colnames:
            injection_catalog["source_type"] = injection_catalog["source_type"].astype(str)

        # Log and return.
        grammar1 = "source" if len(injection_catalog) == 1 else "sources"
        grammar2 = "trixel" if len(injection_catalogs) == 1 else "trixels"
        self.log.info(
            "Retrieved %d injection %s from %d HTM %s.",
            len(injection_catalog),
            grammar1,
            len(injection_catalogs),
            grammar2,
        )
        return injection_catalog

    def _standardize_columns(self, injection_catalog, column_mapping, pixel_scale):
        """Standardize injection catalog column names and units.

        Use config variables to standardize the expected columns and column
        names in the input injection catalog. This method replaces all core
        column names in the config with hard-coded internal names.

        Only a core group of column names are standardized; additional column
        names will not be modified. If certain parameters are needed (i.e.,
        by GalSim), these columns must be given exactly as required in the
        appropriate units. Refer to the configuration documentation for more
        details.

        Parameters
        ----------
        injection_catalog : `astropy.table.Table`
            A catalog of sources to be injected.
        column_mapping : `dict` [`str`, `tuple` [`str`, `astropy.units.Unit`]]
            A dictionary mapping standard column names to the configured column
            names and units.
        pixel_scale : `float`
            Pixel scale of the exposure in arcseconds per pixel.

        Returns
        -------
        injection_catalog : `astropy.table.Table`
            The standardized catalog of sources to be injected.
        """
        self.config = cast(BaseInjectConfig, self.config)

        pixel_scale_equivalency = units.pixel_scale(
            Quantity(pixel_scale, units.arcsec / units.pix)  # type: ignore
        )
        for standard_col, (configured_col, unit) in column_mapping.items():
            # Rename columns if necessary.
            if configured_col in injection_catalog.colnames:
                injection_catalog.rename_column(configured_col, standard_col)
            # Attempt to convert to our desired units, then remove units.
            if standard_col in injection_catalog.columns and unit:
                try:
                    injection_catalog[standard_col] = (
                        injection_catalog[standard_col].to(unit, pixel_scale_equivalency).value
                    )
                except UnitConversionError:
                    pass
        return Table(injection_catalog)

    def _clean_sources(self, injection_catalog, input_exposure):
        """Clean the injection catalog of sources which are not injectable.

        This method will remove sources which are not injectable for a variety
        of reasons, namely: sources which fall outside the padded exposure
        bounding box; sources not selected by virtue of their visit; and
        sources not selected by virtue of their evaluated selection criteria.

        If the input injection catalog contains x/y inputs but does not contain
        RA/Dec inputs, WCS information will be used to generate RA/Dec sky
        coordinate information and appended to the injection catalog.

        Parameters
        ----------
        injection_catalog : `astropy.table.Table`
            The catalog of sources to be injected.
        input_exposure : `lsst.afw.image.exposure.ExposureF`
            The exposure to inject sources into.

        Returns
        -------
        injection_catalog : `astropy.table.Table`
            Updated injection catalog containing *x* and *y* pixel coordinates,
            and cleaned to only include injection sources which fall within the
            bounding box of the input exposure dilated by *trim_padding*.
        """
        self.config = cast(BaseInjectConfig, self.config)

        sources_to_keep = np.ones(len(injection_catalog), dtype=bool)

        # Determine centroids and remove sources outside the padded bbox.
        wcs = input_exposure.getWcs()
        has_sky = {"ra", "dec"} <= set(injection_catalog.columns)
        has_pixel = {"x", "y"} <= set(injection_catalog.columns)
        # Input catalog must contain either RA/Dec OR x/y.
        # If only x/y given, RA/Dec will be calculated.
        if not has_sky and has_pixel:
            begin_x, begin_y = input_exposure.getBBox().getBegin()
            ras, decs = wcs.pixelToSkyArray(
                begin_x + injection_catalog["x"].astype(float),
                begin_y + injection_catalog["y"].astype(float),
                degrees=True,
            )
            injection_catalog["ra"] = ras
            injection_catalog["dec"] = decs
            injection_catalog["x"] += begin_x
            injection_catalog["y"] += begin_y
            has_sky = True
        elif not has_sky and not has_pixel:
            self.log.warning("No spatial coordinates found in injection catalog; cannot inject any sources!")
        if has_sky:
            bbox = input_exposure.getBBox()
            if self.config.trim_padding:
                bbox.grow(int(self.config.trim_padding))
            is_contained = bbox_contains_sky_coords(
                bbox, wcs, injection_catalog["ra"] * units.deg, injection_catalog["dec"] * units.deg
            )
            sources_to_keep &= is_contained
            if (num_not_contained := np.sum(~is_contained)) > 0:
                grammar = ("source", "a centroid") if num_not_contained == 1 else ("sources", "centroids")
                self.log.info(
                    "Removing %d injection %s with %s outside the padded image bounding box.",
                    num_not_contained,
                    grammar[0],
                    grammar[1],
                )

        # Remove sources by boolean selection flag.
        if self.config.selection:
            visit = input_exposure.getInfo().getVisitInfo().getId()
            selected = eval(self.config.selection.format(visit=visit))
            sources_to_keep &= selected
            if (num_not_selected := np.sum(~selected)) >= 0:
                grammar = ["source", "was"] if num_not_selected == 1 else ["sources", "were"]
                self.log.warning(
                    "Removing %d injection %s that %s not selected.",
                    num_not_selected,
                    grammar[0],
                    grammar[1],
                )

        # Print final cleaning report and return.
        num_cleaned_total = np.sum(~sources_to_keep)
        grammar = "source" if len(sources_to_keep) == 1 else "sources"
        self.log.info(
            "Catalog cleaning removed %d of %d %s; %d remaining for catalog flag checking.",
            num_cleaned_total,
            len(sources_to_keep),
            grammar,
            np.sum(sources_to_keep),
        )
        injection_catalog = injection_catalog[sources_to_keep]
        return injection_catalog

    def _flag_sources(self, injection_catalog, input_exposure, binary_flags):
        """Flag sources in the injection catalog prior to source injection.

        This method will flag sources in the injection catalog prior to
        injection. Checks will be made on a number of parameters, including
        magnitude, source type and Sérsic index (where relevant).

        Legacy profile types will be renamed to their standardized GalSim
        equivalents; any source profile types that are not GalSim classes will
        be flagged.

        Note: Unlike the cleaning method, no sources are actually removed here.
        Instead, a binary flag is set in the *injection_flag* column for each
        source. Only unflagged sources will be generated for source injection.

        Parameters
        ----------
        injection_catalog : `astropy.table.Table`
            Catalog of sources to be injected.
        input_exposure : `lsst.afw.image.exposure.ExposureF`
            The exposure sources will be injected into.
        binary_flags : `dict` [`str`, `int`]
            Dictionary of binary flags to be used in the injection_flag column.

        Returns
        -------
        injection_catalog : `astropy.table.Table`
            The cleaned catalog of sources to be injected.
        """
        self.config = cast(BaseInjectConfig, self.config)

        # Flag erroneous magnitude values (missing mag data or NaN mag values).
        if "mag" not in injection_catalog.columns:
            # Check injection_catalog has a mag column.
            self.log.warning("No magnitude data found in injection catalog; cannot inject any sources!")
            injection_catalog["injection_flag"] += 2 ** binary_flags["MAG_BAD"]
        else:
            # Check that all input mag values are finite.
            mag_array = np.isfinite(ma.array(injection_catalog["mag"]))
            bad_mag = ~(mag_array.data * ~mag_array.mask)
            if (num_bad_mag := np.sum(bad_mag)) > 0:
                grammar = "source" if num_bad_mag == 1 else "sources"
                self.log.warning(
                    "Flagging %d injection %s that do not have a finite magnitude.", num_bad_mag, grammar
                )
                injection_catalog["injection_flag"][bad_mag] += 2 ** binary_flags["MAG_BAD"]

        # Replace legacy source types with standardized profile names.
        injection_catalog["source_type"] = injection_catalog["source_type"].astype("O")
        replace_dict = {"Star": "DeltaFunction"}
        for legacy_type, standard_type in replace_dict.items():
            legacy_matches = injection_catalog["source_type"] == legacy_type
            if np.any(legacy_matches):
                injection_catalog["source_type"][legacy_matches] = standard_type
        injection_catalog["source_type"] = injection_catalog["source_type"].astype(str)

        # Flag source types not supported by GalSim.
        input_source_types = set(injection_catalog["source_type"])
        allowed_source_types = [
            "Gaussian",
            "Box",
            "TopHat",
            "DeltaFunction",
            "Airy",
            "Moffat",
            "Kolmogorov",
            "VonKarman",
            "Exponential",
            "DeVaucouleurs",
            "Sersic",
            "InclinedExponential",
            "InclinedSersic",
            "Spergel",
            "RandomKnots",
            "Trail",
            "Stamp",
        ]
        for input_source_type in input_source_types:
            if input_source_type not in allowed_source_types:
                unknown_source_types = injection_catalog["source_type"] == input_source_type
                grammar = "source" if np.sum(unknown_source_types) == 1 else "sources"
                self.log.warning(
                    "Flagging %d injection %s with an unsupported source type: %s.",
                    np.sum(unknown_source_types),
                    grammar,
                    input_source_type,
                )
                injection_catalog["injection_flag"][unknown_source_types] += 2 ** binary_flags["TYPE_UNKNOWN"]

        # Flag extreme Sersic index sources.
        if "n" in injection_catalog.columns:
            min_n = galsim.Sersic._minimum_n
            max_n = galsim.Sersic._maximum_n
            n_vals = injection_catalog["n"]
            extreme_sersics = (n_vals <= min_n) | (n_vals >= max_n)
            if (num_extreme_sersics := np.sum(extreme_sersics)) > 0:
                grammar = "source" if num_extreme_sersics == 1 else "sources"
                self.log.warning(
                    "Flagging %d injection %s with a Sersic index outside the range %.1f <= n <= %.1f.",
                    num_extreme_sersics,
                    grammar,
                    min_n,
                    max_n,
                )
                injection_catalog["injection_flag"][extreme_sersics] += 2 ** binary_flags["SERSIC_EXTREME"]

        # Print final cleaning report.
        num_flagged_total = np.sum(injection_catalog["injection_flag"] != 0)
        grammar = "source" if len(injection_catalog) == 1 else "sources"
        self.log.info(
            "Catalog flag checking applied flags to %d of %d %s; %d remaining for source generation.",
            num_flagged_total,
            len(injection_catalog),
            grammar,
            np.sum(injection_catalog["injection_flag"] == 0),
        )
        return injection_catalog
