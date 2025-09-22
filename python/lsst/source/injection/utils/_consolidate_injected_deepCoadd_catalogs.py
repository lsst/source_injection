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

__all__ = [
    "ConsolidateInjectedCatalogsConnections",
    "ConsolidateInjectedCatalogsConfig",
    "ConsolidateInjectedCatalogsTask",
]

from collections import defaultdict

import astropy.table
import astropy.units as u
import numpy as np
from astropy.table import Table, join, vstack
from astropy.table.column import MaskedColumn
from smatch.matcher import Matcher  # type: ignore [import-not-found]

from lsst.daf.butler import DatasetProvenance
from lsst.geom import Box2D, SpherePoint, degrees
from lsst.pex.config import Field, ListField
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base.connections import InputQuantizedConnection
from lsst.pipe.base.connectionTypes import Input, Output
from lsst.skymap import BaseSkyMap


class ConsolidateInjectedCatalogsConnections(
    PipelineTaskConnections,
    dimensions=("instrument", "skymap", "tract"),
    defaultTemplates={
        "injected_prefix": "injected_",
    },
):
    """Base connections for source injection tasks."""

    input_catalogs = Input(
        doc="Per-patch and per-band injected catalogs to draw inputs from.",
        name="{injected_prefix}deepCoadd_catalog",
        dimensions=("skymap", "tract", "patch", "band"),
        storageClass="ArrowAstropy",
        minimum=1,
        multiple=True,
    )
    skyMap = Input(
        doc="Input definition of geometry/bbox and projection/wcs for warped exposures",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    output_catalog = Output(
        doc="Per-tract multiband catalog of injected sources.",
        name="{injected_prefix}deepCoadd_catalog_tract",
        storageClass="ArrowAstropy",
        dimensions=("skymap", "tract"),
    )


def _get_catalogs(
    inputs: dict,
    input_refs: InputQuantizedConnection,
    skymap: BaseSkyMap,
    col_ra: str = "ra",
    col_dec: str = "dec",
    include_outer: bool = True,
) -> tuple[dict, int]:
    """Organize input catalogs into a dictionary with photometry band
    keys.

    Parameters
    ----------
    inputs: `dict`
        A dictionary containing the input datasets.
    input_refs: `lsst.pipe.base.connections.InputQuantizedConnection`
        The input dataset references used by the butler.
    skymap : `lsst.skymap.BaseSkyMap`
        Sky tessellation object
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).
    include_outer: `bool`
        Whether to include objects injected into the outer (not inner) region
        of a patch.

    Returns
    -------
    `tuple[dict, int]`
        contains :
            catalog_dict: `dict`
                A dictionary with photometric bands for keys and astropy
                tables for items.
            tract: `int`
                The tract covering the catalogs in catalog_dict.
    """
    catalog_dict: dict[str, dict[str, astropy.table.Table]] = {}
    tracts = set()
    for ref, catalog in zip(input_refs.input_catalogs, inputs["input_catalogs"]):
        band = ref.dataId.band.name
        if band not in catalog_dict:
            catalog_dict[band] = {}
        tract = ref.dataId.tract.id
        tractInfo = skymap[tract]
        # Load the patch number to check for patch overlap duplicates.
        patch = ref.dataId.patch.id
        if not include_outer:
            is_inner = getPatchInner(
                tractInfo[patch],
                catalog[col_ra],
                catalog[col_dec],
            )
            catalog = catalog[is_inner]
        catalog["patch"] = patch
        # Strip provenance from catalogs before merging to avoid the
        # provenance headers triggering warnings in the astropy naive
        # metadata merge tool.
        DatasetProvenance.strip_provenance_from_flat_dict(catalog.meta)
        catalog_dict[band][patch] = catalog
        tracts.add(tract)
        # Check that only catalogs covered by a single tract are loaded.
        if len(tracts) > 1:
            raise RuntimeError(f"Got tract={tract} with {tracts=}; there should only be one tract")
    # Stack the per-band catalogs.
    for band, catalog_patches in catalog_dict.items():
        catalog_dict[band] = vstack([item[1] for item in sorted(catalog_patches.items())])
    return catalog_dict, list(tracts)[0]


def _get_patches(
    catalog_dict: dict[str, astropy.table.Table],
    tractInfo,
    col_ra: str = "ra",
    col_dec: str = "dec",
    index=None,
):
    """Create a patch column and assign each row a patch number.

    Parameters
    ----------
    catalog_dict: `dict`
        A dictionary with photometric bands for keys and astropy tables for
        items.
    tractInfo: `lsst.skymap.tractInfo.ExplicitTractInfo`
        Information for a tract specified explicitly.
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).
    """
    for catalog in catalog_dict.values():
        if "patch" not in catalog.colnames:
            patches = np.empty(len(catalog), dtype=int)
            catalog.add_column(patches, name="patch", index=index)

        patches = catalog["patch"]
        for idx, row in enumerate(catalog):
            coord = SpherePoint(row[col_ra], row[col_dec], degrees)
            patchInfo = tractInfo.findPatch(coord)
            patches[idx] = int(patchInfo.getSequentialIndex())


def getPatchInner(
    patchInfo,
    ra: np.ndarray,
    dec: np.ndarray,
):
    """Set a flag for each source if it is in the innerBBox of a patch.

    Parameters
    ----------
    patchInfo : `lsst.skymap.PatchInfo`
        Information about a `SkyMap` `Patch`.
    ra: `np.ndarray`
        Right ascension values in degrees.
    dec: `np.ndarray`
        Declination values in degrees.

    Returns
    -------
    isPatchInner : array-like of `bool`
        `True` for each source that has a centroid
        in the inner region of a patch.
    """
    # convert the coordinates to pixel positions
    wcs = patchInfo.getWcs()
    x, y = wcs.skyToPixelArray(ra, dec, degrees=True)

    # set inner flags for each source
    innerFloatBBox = Box2D(patchInfo.getInnerBBox())
    isPatchInner = innerFloatBBox.contains(x, y)

    return isPatchInner


def getTractInner(
    tractInfo,
    skyMap,
    ra: np.ndarray,
    dec: np.ndarray,
):
    """Set a flag for each source that the skyMap includes in tractInfo.

    Parameters
    ----------
    tractInfo : `lsst.skymap.TractInfo`
        Tract object
    skyMap : `lsst.skymap.BaseSkyMap`
        Sky tessellation object
    ra: `np.ndarray`
        Right ascension values in degrees.
    dec: `np.ndarray`
        Declination values in degrees.

    Returns
    -------
    isTractInner : array-like of `bool`
        True if the skyMap.findTract method returns
        the same tract as tractInfo.
    """
    tractId = tractInfo.getId()
    isTractInner = np.array(
        [skyMap.findTract(SpherePoint(_ra, _dec, degrees)).getId() == tractId for _ra, _dec in zip(ra, dec)],
        dtype=bool,
    )

    return isTractInner


class ConsolidateInjectedCatalogsConfig(  # type: ignore [call-arg]
    PipelineTaskConfig, pipelineConnections=ConsolidateInjectedCatalogsConnections
):
    """Base configuration for source injection tasks."""

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
        "brightness profiles defined by GalSim. For more information see the Galsim docs at "
        "https://galsim-developers.github.io/GalSim/_build/html/sb.html",
        default="source_type",
    )
    columns_extra = ListField[str](
        doc="Extra columns to be copied from the injection catalog (e.g. for shapes)",
        default=[],
    )
    groupIdKey = Field[str](
        doc="Key for the group id column to merge sources on, if any",
        default=None,
        optional=True,
    )
    injectionKey = Field[str](
        doc="True if the source was successfully injected.",
        default="injection_flag",
    )
    injectionSizeKey = Field[str](
        doc="The size of the drawn injection box",
        default="injection_draw_size",
    )
    isPatchInnerKey = Field[str](
        doc="True if source is in the inner region of a coadd patch.",
        default="injected_isPatchInner",
    )
    isTractInnerKey = Field[str](
        doc="True if source is in the inner region of a coadd tract.",
        default="injected_isTractInner",
    )
    isPrimaryKey = Field[str](
        doc="True if the source was successfully injected and is in both the inner region of a coadd patch "
        "and tract.",
        default="injected_isPrimary",
    )
    remove_patch_overlap_duplicates = Field[bool](
        doc="Optional parameter to remove patch overlap duplicate sources.",
        default=False,
    )
    get_catalogs_from_butler = Field[bool](
        doc="Optional parameter to specify whether or not the input catalogs are loaded with a Butler.",
        default=True,
    )
    pixel_match_radius = Field[float](
        doc="Radius for matching catalogs across different bands.",
        default=0.1,
    )

    def consolidate_deepCoadd(
        self,
        catalog_dict: dict[str, astropy.table.Table],
        skymap: BaseSkyMap,
        tract: int,
        copy_catalogs: bool = False,
    ) -> Table:
        """Consolidate all tables in catalog_dict into one table.

        Parameters
        ----------
        catalog_dict: `dict`
            A dictionary with photometric bands for keys and astropy tables for
            items.
        skymap: `lsst.skymap.BaseSkyMap`
            A base skymap.
        tract: `int`
            The tract where sources have been injected.
        copy_catalogs: `bool`
            Whether to copy the input catalogs; if False, they will be modified
            in-place.

        Returns
        -------
        multiband_catalog: `astropy.table.Table`
            A single table containing all information of the separate
            tables in catalog_dict
        """
        tractInfo = skymap.generateTract(tract)
        # If patch numbers are not loaded via dataIds from the butler, manually
        # load patch numbers from source positions.
        if not self.get_catalogs_from_butler:
            _get_patches(catalog_dict, tractInfo)

        # Convert the pixel match radius to degrees.
        tractWcs = tractInfo.getWcs()
        pixel_scale = tractWcs.getPixelScale()
        match_radius = self.pixel_match_radius * pixel_scale.asDegrees()
        has_groups = bool(self.groupIdKey)
        bands = list(catalog_dict.keys())
        if has_groups or (len(bands) > 1):
            # Match the catalogs across bands.
            output_catalog = self.make_multiband_catalog(
                bands,
                catalog_dict,
                match_radius,
                copy_catalogs=copy_catalogs,
            )
        else:
            output_catalog = catalog_dict[bands[0]]
            output_catalog.rename_column(self.col_mag, f"{bands[0]}_{self.col_mag}")
        # Remove sources outside tract boundaries.
        out_of_tract_bounds = []
        for index, (ra, dec) in enumerate(zip(output_catalog[self.col_ra], output_catalog[self.col_dec])):
            point = SpherePoint(ra * degrees, dec * degrees)
            if not tractInfo.contains(point):
                out_of_tract_bounds.append(index)
        if out_of_tract_bounds:
            output_catalog.remove_rows(out_of_tract_bounds)
        # Assign flags.
        # There may be a pre-existing patch column; however, patches can shift
        # if centroids vary per band, etc.
        # TODO: Need to check if this can cause duplicate or dropped entries
        # for objects near patch boundaries
        _get_patches(
            {"x": output_catalog},
            tractInfo=tractInfo,
            col_ra=self.col_ra,
            col_dec=self.col_dec,
        )
        patches = list(set(output_catalog["patch"]))
        self.setPrimaryFlags(
            catalog=output_catalog,
            skyMap=skymap,
            tractInfo=tractInfo,
            patches=patches,
        )
        # Add a new injected_id column.
        output_catalog.add_column(col=np.arange(len(output_catalog)), name="injected_id")
        # If using a group ID, the draw_size column is
        # not preserved, and the rest are already ordered
        if self.groupIdKey:
            return output_catalog
        # Reorder columns
        mag_cols = [col for col in output_catalog.columns if f"_{self.col_mag}" in col]

        new_order = [
            "injected_id",
            self.col_ra,
            self.col_dec,
            self.col_source_type,
            *mag_cols,
            "patch",
            "injection_id",
            self.injectionSizeKey,
            self.injectionKey,
            "injected_isPatchInner",
            "injected_isTractInner",
            "injected_isPrimary",
        ]
        for column in output_catalog.columns:
            if column not in new_order:
                new_order.append(column)

        return output_catalog[new_order]

    def make_multiband_catalog(
        self,
        bands: list,
        catalog_dict: dict[str, astropy.table.Table],
        match_radius: float,
        copy_catalogs: bool = False,
    ) -> Table:
        """Combine multiple band-specific catalogs into one multiband
        catalog.

        Parameters
        ----------
        bands: `list`
            A list of string photometry bands.
        catalog_dict: `dict`
            A dictionary with photometric bands for keys and astropy
            tables for items.
        match_radius: `float`
            The radius for matching catalogs across bands in arcsec.
        copy_catalogs: `bool`
            Whether to copy the input catalogs; if False, they will be modified
            in-place.

        Returns
        -------
        multiband_catalog: `astropy.table.Table`
            A catalog with sources that have magnitude information across all
            bands.
        """
        col_mag = self.col_mag
        col_ra = self.col_ra
        col_dec = self.col_dec

        if self.groupIdKey:
            n_comps = 0
            n_rows = {}
            for band, catalog in catalog_dict.items():
                groupIds, counts = np.unique(
                    catalog[self.groupIdKey],
                    return_counts=True,
                )
                n_rows[band] = len(groupIds)
                n_comps = max(np.max(counts), n_comps)

            # Maybe validate that this is true? Why set groupId otherwise?
            is_multicomp = n_comps > 1
            prefix_comp = "comp{idx_comp}_" if is_multicomp else ""
            columns_comp = {
                column_in: f"{{band}}_{prefix_comp}{column_in}"
                for column_in in [self.col_source_type, self.injectionKey] + list(self.columns_extra)
            }

            for band, catalog in catalog_dict.items():
                n_rows_b = n_rows[band]
                unit_mag = catalog["mag"].unit or u.ABmag
                counts = defaultdict(int)
                idxs_new = {}

                groupIds = np.full(n_rows_b, 0, dtype=catalog[self.groupIdKey].dtype)
                ra = np.full(n_rows_b, np.nan, dtype=catalog[col_ra].dtype)
                dec = np.full(n_rows_b, np.nan, dtype=catalog[col_dec].dtype)
                injection_flag = np.full(n_rows_b, True, dtype=bool)
                if is_multicomp:
                    flux = np.full(n_rows_b, np.nan, dtype=catalog[col_mag].dtype)
                else:
                    mag = np.full(n_rows_b, np.nan, dtype=catalog[col_mag].dtype)

                values_comp = {}
                if is_multicomp:
                    for idx_comp in range(1, n_comps + 1):
                        prefix_comp = f"comp{idx_comp}_"
                        values_comp[f"{band}_{prefix_comp}flux"] = np.ma.masked_array(
                            np.full(n_rows_b, np.nan, dtype=catalog[col_mag].dtype),
                            mask=np.ones(n_rows_b, dtype=bool),
                        )
                for idx_comp in range(1, n_comps + 1):
                    for column_comp_in, column_comp_out in columns_comp.items():
                        values_col = np.full(n_rows_b, np.nan, dtype=catalog[column_comp_in].dtype)
                        column_out = column_comp_out.format(band=band, idx_comp=idx_comp)
                        values_comp[column_out] = (
                            np.ma.masked_array(values_col, mask=np.ones(n_rows_b, dtype=bool))
                            if is_multicomp
                            else values_col
                        )

                idx_new = 0
                for row in catalog:
                    groupId = row[self.groupIdKey]
                    idx_comp = counts[groupId] + 1
                    counts[groupId] = idx_comp
                    prefix_comp = f"comp{idx_comp}_" if is_multicomp else ""

                    injected = row[self.injectionKey] == False  # noqa: E712
                    if idx_comp == 1:
                        ra[idx_new] = row[col_ra]
                        dec[idx_new] = row[col_dec]
                        groupIds[idx_new] = groupId
                        idxs_new[groupId] = idx_new
                        if is_multicomp:
                            flux_comp = (row[col_mag] * unit_mag).to(u.nJy).value
                            flux[idx_new] = flux_comp * injected
                        else:
                            mag[idx_new] = row[col_mag]
                        idx_old = idx_new
                        idx_new += 1
                    else:
                        idx_old = idxs_new[groupId]
                        flux_cumul = flux[idx_old]
                        flux_comp = (row[col_mag] * unit_mag).to(u.nJy).value
                        # TODO: Excluding not-injected components means
                        # objects with no injected components will have the
                        # ra,dec of their first component, not the mean of
                        # all excluded components.
                        if flux_comp * injected > 0:
                            flux_new = flux_cumul + flux_comp
                            flux[idx_old] = flux_new
                            if ra[idx_old] != row[col_ra]:
                                # Take a weighted mean
                                # TODO: deal with periodicity
                                ra[idx_old] = (flux_cumul * ra[idx_old] + flux_comp * row[col_ra]) / flux_new
                            if dec[idx_old] != row[col_dec]:
                                dec[idx_old] = (
                                    flux_cumul * dec[idx_old] + flux_comp * row[col_dec]
                                ) / flux_new

                    # One component is sufficient to say it was injected
                    injection_flag[idx_old] &= not injected

                    if is_multicomp:
                        column = values_comp[f"{band}_{prefix_comp}flux"]
                        column.mask[idx_old] = False
                        column[idx_old] = flux_comp

                    for column_comp in columns_comp:
                        column = f"{band}_{prefix_comp}{column_comp}"
                        values_comp[column].mask[idx_old] = False
                        values_comp[column][idx_old] = row[column_comp]

                if is_multicomp:
                    mag = (flux * u.nJy).to(u.ABmag).value

                columns = {
                    self.groupIdKey: groupIds,
                    f"{band}_{self.injectionKey}": injection_flag,
                    f"{band}_{col_ra}": ra,
                    f"{band}_{col_dec}": dec,
                    f"{band}_{col_mag}": mag,
                }
                units = {
                    self.groupIdKey: None,
                    f"{band}_{self.injectionKey}": None,
                    f"{band}_{col_ra}": catalog[col_ra].unit,
                    f"{band}_{col_dec}": catalog[col_dec].unit,
                    f"{band}_{col_mag}": catalog[col_mag].unit,
                }

                for idx_comp in range(1, n_comps + 1):
                    for column_comp_in, column_comp_out in columns_comp.items():
                        name_column = column_comp_out.format(band=band, idx_comp=idx_comp)
                        columns[name_column] = values_comp[name_column]
                        units[name_column] = catalog[column_comp_in].unit

                catalog_new = astropy.table.Table(columns, units=units)
                catalog_dict[band] = catalog_new

            copy_catalogs = False

        # Load the first catalog then loop to add info for the other bands.
        multiband_catalog = catalog_dict[bands[0]]
        if copy_catalogs:
            multiband_catalog = multiband_catalog.copy()
        if not self.groupIdKey:
            multiband_catalog.add_column(col=multiband_catalog[col_mag], name=f"{bands[0]}_{col_mag}")
            multiband_catalog.remove_column(col_mag)
        else:
            coords_same = True

        for band in bands[1:]:
            if not self.groupIdKey:
                # Make a column for the new band.
                multiband_catalog.add_column([np.nan] * len(multiband_catalog), name=f"{band}_{col_mag}")
            # Match the input catalog for this band to the existing
            # multiband catalog.
            catalog_next_band = catalog_dict[band]
            if copy_catalogs:
                catalog_next_band = catalog_next_band.copy()
            if not self.groupIdKey:
                catalog_next_band.rename_column(col_mag, f"{band}_{col_mag}")
            if match_radius >= 0:
                with Matcher(multiband_catalog[col_ra], multiband_catalog[col_dec]) as m:
                    idx, multiband_match_inds, next_band_match_inds, dists = m.query_radius(
                        catalog_next_band[col_ra],
                        catalog_next_band[col_dec],
                        match_radius,
                        return_indices=True,
                    )
            else:
                if self.groupIdKey:
                    multiband_catalog = join(
                        multiband_catalog,
                        catalog_next_band,
                        keys=self.groupIdKey,
                        join_type="outer",
                    )
                    if coords_same:
                        for coord in (col_ra, col_dec):
                            coords = multiband_catalog[f"{band}_{coord}"]
                            coords_good = coords[np.isfinite(coords)]
                            coords_good_first = multiband_catalog[f"{bands[0]}_{coord}"]
                            coords_good_first = coords_good_first[np.isfinite(coords_good_first)]
                            if (len(coords_good) != len(coords_good_first)) or (
                                np.any(coords_good != coords_good_first)
                            ):
                                coords_same = False

                multiband_match_inds, next_band_match_inds = [], []
            # If there are matches...
            if len(multiband_match_inds) > 0 and len(next_band_match_inds) > 0:
                # ...choose the coordinates in the brightest band.
                for i, j in zip(multiband_match_inds, next_band_match_inds):
                    mags = []
                    for col in multiband_catalog.colnames:
                        if f"_{col_mag}" in col:
                            mags.append((col, multiband_catalog[i][col]))
                    bright_mag = min([x[1] for x in mags])
                    if catalog_next_band[f"{band}_{col_mag}"][j] < bright_mag:
                        multiband_catalog[col_ra][i] = catalog_next_band[col_ra][j]
                        multiband_catalog[col_dec][i] = catalog_next_band[col_dec][j]
                # TODO: Once multicomponent object support is added, make some
                # logic to pick the correct source_type.
                # Fill the new mag value.
                multiband_catalog[f"{band}_{col_mag}"][multiband_match_inds] = catalog_next_band[
                    f"{band}_{col_mag}"
                ][next_band_match_inds]
                # Add rows for all the sources without matches.
                not_next_band_match_inds = np.full(len(catalog_next_band), True, dtype=bool)
                not_next_band_match_inds[next_band_match_inds] = False
                multiband_catalog = vstack([multiband_catalog, catalog_next_band[not_next_band_match_inds]])
            # Otherwise just stack the tables.
            elif not self.groupIdKey:
                multiband_catalog = vstack([multiband_catalog, catalog_next_band])

        if self.groupIdKey:
            if coords_same:
                for coord in (col_ra, col_dec):
                    for band in bands[1:]:
                        del multiband_catalog[f"{band}_{coord}"]
                    multiband_catalog.rename_column(f"{bands[0]}_{coord}", coord)
            else:
                fluxes = np.array(
                    [(multiband_catalog[f"{band}_{col_mag}"] * u.ABmag).to(u.nJy).value for band in bands]
                )
                # TODO: Test this better and deal with periodicity in RA
                for coord in (col_dec, col_ra):
                    coords = np.array([multiband_catalog[f"{band}_{coord}"] for band in bands])
                    coords = np.nanmean(fluxes * coords, axis=0) / np.nansum(fluxes, axis=0)
                    multiband_catalog.add_column(coords, index=1, name=coord)
            multiband_catalog.add_column(
                np.all(
                    np.array([multiband_catalog[f"{band}_{self.injectionKey}"] for band in bands]) == 1,
                    axis=0,
                ),
                index=1,
                name=self.injectionKey,
            )
        else:
            # Fill in per-band injection flag columns
            if not copy_catalogs:
                multiband_catalog[self.injectionKey][:] = catalog_dict[bands[0]][self.injectionKey][:]
            multiband_catalog[f"{bands[0]}_{self.injectionKey}"] = multiband_catalog[self.injectionKey]
            for band in bands[1:]:
                injected_band = catalog_dict[band][self.injectionKey]
                multiband_catalog[self.injectionKey] &= injected_band
                multiband_catalog[f"{band}_{self.injectionKey}"] = injected_band

        # Fill any automatically masked values with NaNs if possible (float)
        # Otherwise, use the dtype's minimum value (for int, bool, etc.)
        if multiband_catalog.has_masked_columns:
            for colname in multiband_catalog.columns:
                column = multiband_catalog[colname]
                if isinstance(column, MaskedColumn):
                    # Set the underlying values in-place
                    column._data[column.mask] = (
                        np.nan if np.issubdtype(column.dtype, float) else np.ma.maximum_fill_value(column)
                    )
        return multiband_catalog

    def setPrimaryFlags(
        self,
        catalog,
        skyMap,
        tractInfo,
        patches: list,
    ):
        """Set isPrimary and related flags on sources.

        For co-added imaging, the `isPrimary` flag returns True when an object
        has no children, is in the inner region of a coadd patch, is in the
        inner region of a coadd tract, and is not detected in a pseudo-filter
        (e.g., a sky_object).
        For single frame imaging, the isPrimary flag returns True when a
        source has no children and is not a sky source.

        Parameters
        ----------
        catalog: `astropy.table.Table`
            A catalog of sources.
            Writes is-patch-inner, is-tract-inner, and is-primary flags.
        skyMap : `lsst.skymap.BaseSkyMap`
            Sky tessellation object
        tractInfo : `lsst.skymap.TractInfo`
            Tract object
        patches : `list`
            List of coadd patches
        """
        # Mark whether sources are contained within the inner regions of the
        # given tract/patch.
        isPatchInner = np.full(len(catalog), 0, dtype=bool)
        ra, dec = catalog[self.col_ra].data, catalog[self.col_dec].data
        for patch in patches:
            patchMask = catalog["patch"] == patch
            patchInfo = tractInfo.getPatchInfo(patch)
            isPatchInner[patchMask] = getPatchInner(patchInfo, ra[patchMask], dec[patchMask])
        isTractInner = getTractInner(tractInfo, skyMap, ra, dec)
        isPrimary = isTractInner & isPatchInner & (catalog[self.injectionKey] == 0)

        catalog[self.isPatchInnerKey] = isPatchInner
        catalog[self.isTractInnerKey] = isTractInner
        catalog[self.isPrimaryKey] = isPrimary


class ConsolidateInjectedCatalogsTask(PipelineTask):
    """Class for combining all tables in a collection of input catalogs
    into one table.
    """

    _DefaultName = "consolidateInjectedCatalogsTask"
    ConfigClass = ConsolidateInjectedCatalogsConfig

    def runQuantum(self, butlerQC, input_refs, output_refs):
        inputs = butlerQC.get(input_refs)
        skymap = inputs["skyMap"]
        catalog_dict, tract = _get_catalogs(
            inputs,
            input_refs,
            skymap,
            col_ra=self.config.col_ra,
            col_dec=self.config.col_dec,
            include_outer=not bool(self.config.groupIdKey),
        )
        outputs = self.run(catalog_dict, skymap, tract)
        butlerQC.put(outputs, output_refs)

    def run(
        self,
        catalog_dict: dict,
        skymap: BaseSkyMap,
        tract: int,
    ) -> Table:
        """Consolidate all tables in catalog_dict into one table.

        catalog_dict: `dict`
            A dictionary with photometric bands for keys and astropy tables for
            items.
        skymap: `lsst.skymap.BaseSkyMap`
            A base skymap.
        tract: `int`
            The tract where sources have been injected.

        Returns
        -------
        output_struct : `lsst.pipe.base.Struct`
            contains :
                multiband_catalog: `astropy.table.Table`
                    A single table containing all information of the separate
                    tables in catalog_dict
        """
        output_catalog = self.config.consolidate_deepCoadd(
            catalog_dict,
            skymap,
            tract,
        )
        output_struct = Struct(output_catalog=output_catalog)
        return output_struct
