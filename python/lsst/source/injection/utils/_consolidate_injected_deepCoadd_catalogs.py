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
    "consolidate_injected_deepCoadd_catalogs",
]

import numpy as np
from astropy.table import Table, vstack
from astropy.table.column import MaskedColumn
from lsst.geom import Box2D, SpherePoint, degrees
from lsst.pex.config import Field
from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections, Struct
from lsst.pipe.base.connections import InputQuantizedConnection
from lsst.pipe.base.connectionTypes import Input, Output, PrerequisiteInput
from lsst.skymap import BaseSkyMap
from smatch.matcher import Matcher  # type: ignore [import-not-found]


class ConsolidateInjectedCatalogsConnections(
    PipelineTaskConnections,
    dimensions=("instrument", "skymap", "tract"),
    defaultTemplates={
        "injected_prefix": "injected_",
    },
):
    """Base connections for source injection tasks."""

    input_catalogs = PrerequisiteInput(
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
    injectionKey = Field[str](
        doc="True if the source was successfully injected.",
        default="injection_flag",
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


def _get_catalogs(
    inputs: dict,
    input_refs: InputQuantizedConnection,
) -> tuple[dict, int]:
    """Organize input catalogs into a dictionary with photometry band
    keys.

    Parameters
    ----------
    inputs: `dict`
        A dictionary containing the input datasets.
    input_refs: `lsst.pipe.base.connections.InputQuantizedConnection`
        The input dataset references used by the butler.

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
    catalog_dict: dict[str, list] = {}
    tracts = set()
    for ref, catalog in zip(input_refs.input_catalogs, inputs["input_catalogs"]):
        band = ref.dataId.band.name
        if band not in catalog_dict:
            catalog_dict[band] = []
        # Load the patch number to check for patch overlap duplicates.
        catalog["patch"] = ref.dataId.patch.id
        catalog_dict[band].append(catalog)
        tracts.add(ref.dataId.tract.id)
    # Stack the per-band catalogs.
    for band, catalog_list in catalog_dict.items():
        catalog_dict[band] = vstack(catalog_list)
    # Check that only catalogs covered by a single tract are loaded.
    if len(tracts) != 1:
        raise RuntimeError(f"len({tracts=}) != 1")
    return (catalog_dict, list(tracts)[0])


def _get_patches(
    catalog_dict: dict,
    tractInfo,
    col_ra: str = "ra",
    col_dec: str = "dec",
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
        if "patch" not in list(set(catalog.columns)):
            catalog.add_column(col=-9999, name="patch")
        for row in catalog:
            coord = SpherePoint(row[col_ra], row[col_dec], degrees)
            patchInfo = tractInfo.findPatch(coord)
            row["patch"] = int(patchInfo.getSequentialIndex())


def getPatchInner(
    catalog: Table,
    patchInfo,
    col_ra: str = "ra",
    col_dec: str = "dec",
):
    """Set a flag for each source if it is in the innerBBox of a patch.

    Parameters
    ----------
    catalog: `astropy.table.Table`
        A catalog of sources.
    patchInfo : `lsst.skymap.PatchInfo`
        Information about a `SkyMap` `Patch`.
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).

    Returns
    -------
    isPatchInner : array-like of `bool`
        `True` for each source that has a centroid
        in the inner region of a patch.
    """
    # Extract positions for all the sources.
    ra = catalog[col_ra]
    dec = catalog[col_dec]

    # convert the coordinates to pixel positions
    wcs = patchInfo.getWcs()
    x, y = wcs.skyToPixelArray(ra, dec, degrees=True)

    # set inner flags for each source
    innerFloatBBox = Box2D(patchInfo.getInnerBBox())
    isPatchInner = innerFloatBBox.contains(x, y)

    return isPatchInner


def getTractInner(
    catalog,
    tractInfo,
    skyMap,
    col_ra: str = "ra",
    col_dec: str = "dec",
):
    """Set a flag for each source that the skyMap includes in tractInfo.

    Parameters
    ----------
    catalog: `astropy.table.Table`
        A catalog of sources.
    tractInfo : `lsst.skymap.TractInfo`
        Tract object
    skyMap : `lsst.skymap.BaseSkyMap`
        Sky tessellation object
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).

    Returns
    -------
    isTractInner : array-like of `bool`
        True if the skyMap.findTract method returns
        the same tract as tractInfo.
    """
    ras, decs = catalog[col_ra].value, catalog[col_dec].value
    skyCoords = [SpherePoint(ra, dec, degrees) for (ra, dec) in list(zip(ras, decs))]
    tractId = tractInfo.getId()
    isTractInner = np.array([skyMap.findTract(coord).getId() == tractId for coord in skyCoords])
    return isTractInner


def setPrimaryFlags(
    catalog,
    skyMap,
    tractInfo,
    patches: list,
    col_ra: str = "ra",
    col_dec: str = "dec",
    isPatchInnerKey: str = "injected_isPatchInner",
    isTractInnerKey: str = "injected_isTractInner",
    isPrimaryKey: str = "injected_isPrimary",
    injectionKey: str = "injection_flag",
):
    """Set isPrimary and related flags on sources.

    For coadded imaging, the `isPrimary` flag returns True when an object
    has no children, is in the inner region of a coadd patch, is in the
    inner region of a coadd trach, and is not detected in a pseudo-filter
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
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).
    isPatchInnerKey: `str`
        Column name for the isPatchInner flag.
    isTractInnerKey: `str`
        Column name for the isTractInner flag.
    isPrimaryKey: `str`
        Column name for the isPrimary flag.
    injectionKey: `str`
        Column name for the injection flag.
    """
    # Mark whether sources are contained within the inner regions of the
    # given tract/patch.
    isPatchInner = np.array([False] * len(catalog))
    for patch in patches:
        patchMask = catalog["patch"] == patch
        patchInfo = tractInfo.getPatchInfo(patch)
        isPatchInner[patchMask] = getPatchInner(catalog[patchMask], patchInfo, col_ra, col_dec)
    isTractInner = getTractInner(catalog, tractInfo, skyMap)
    isPrimary = isTractInner & isPatchInner & (catalog[injectionKey] == 0)

    catalog[isPatchInnerKey] = isPatchInner
    catalog[isTractInnerKey] = isTractInner
    catalog[isPrimaryKey] = isPrimary


def _make_multiband_catalog(
    bands: list,
    catalog_dict: dict,
    match_radius: float,
    col_ra: str = "ra",
    col_dec: str = "dec",
    col_mag: str = "mag",
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
        The radius for matching catalogs across bands.
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).
    col_mag: `str`
        Column name for magnitude.

    Returns
    -------
    multiband_catalog: `astropy.table.Table`
        A catalog with sources that have magnitude information across all
        bands.
    """
    # Load the first catalog then loop to add info for the other bands.
    multiband_catalog = catalog_dict[bands[0]].copy()
    multiband_catalog.add_column(col=multiband_catalog[col_mag], name=f"{bands[0]}_{col_mag}")
    multiband_catalog.remove_column(col_mag)
    for band in bands[1:]:
        # Make a column for the new band.
        multiband_catalog.add_column([np.nan] * len(multiband_catalog), name=f"{band}_{col_mag}")
        # Match the input catalog for this band to the existing
        # multiband catalog.
        catalog_next_band = catalog_dict[band].copy()
        catalog_next_band.rename_column(col_mag, f"{band}_{col_mag}")
        with Matcher(multiband_catalog[col_ra], multiband_catalog[col_dec]) as m:
            idx, multiband_match_inds, next_band_match_inds, dists = m.query_radius(
                catalog_next_band[col_ra],
                catalog_next_band[col_dec],
                match_radius,
                return_indices=True,
            )
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
        else:
            multiband_catalog = vstack([multiband_catalog, catalog_next_band])
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


def consolidate_injected_deepCoadd_catalogs(
    catalog_dict: dict,
    skymap: BaseSkyMap,
    tract: int,
    pixel_match_radius: float = 0.1,
    get_catalogs_from_butler: bool = True,
    col_ra: str = "ra",
    col_dec: str = "dec",
    col_mag: str = "mag",
    isPatchInnerKey="injected_isPatchInner",
    isTractInnerKey="injected_isTractInner",
    isPrimaryKey="injected_isPrimary",
    injectionKey="injection_flag",
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
    pixel_match_radius: `float`
        Match radius in pixels to use for self-matching catalogs across
        different bands.
    get_catalogs_from_butler: `bool`
        Optional parameter to specify whether or not the input catalogs are
        loaded with a Butler.
    col_ra: `str`
        Column name for right ascension (in degrees).
    col_dec: `str`
        Column name for declination (in degrees).
    col_mag: `str`
        Column name for magnitude.
    isPatchInnerKey: `str`
        Column name for the isPatchInner flag.
    isTractInnerKey: `str`
        Column name for the isTractInner flag.
    isPrimaryKey: `str`
        Column name for the isPrimary flag.
    injectionKey: `str`
        Column name for the injection flag.

    Returns
    -------
    multiband_catalog: `astropy.table.Table`
        A single table containing all information of the separate
        tables in catalog_dict
    """
    tractInfo = skymap.generateTract(tract)
    # If patch numbers are not loaded via dataIds from the butler, manually
    # load patch numbers from source positions.
    if not get_catalogs_from_butler:
        _get_patches(catalog_dict, tractInfo)

    # Convert the pixel match radius to degrees.
    tractWcs = tractInfo.getWcs()
    pixel_scale = tractWcs.getPixelScale()
    match_radius = pixel_match_radius * pixel_scale.asDegrees()
    bands = list(catalog_dict.keys())
    if len(bands) > 1:
        # Match the catalogs across bands.
        output_catalog = _make_multiband_catalog(bands, catalog_dict, match_radius)
    else:
        output_catalog = catalog_dict[bands[0]]
        output_catalog.rename_column(col_mag, f"{bands[0]}_{col_mag}")
    # Remove sources outside tract boundaries.
    out_of_tract_bounds = []
    for index, (ra, dec) in enumerate(list(zip(output_catalog[col_ra], output_catalog[col_dec]))):
        point = SpherePoint(ra * degrees, dec * degrees)
        if not tractInfo.contains(point):
            out_of_tract_bounds.append(index)
    output_catalog.remove_rows(out_of_tract_bounds)
    # Assign flags.
    patches = list(set(output_catalog["patch"]))
    setPrimaryFlags(
        catalog=output_catalog,
        skyMap=skymap,
        tractInfo=tractInfo,
        patches=patches,
        col_ra=col_ra,
        col_dec=col_dec,
        isPatchInnerKey=isPatchInnerKey,
        isTractInnerKey=isTractInnerKey,
        isPrimaryKey=isPrimaryKey,
        injectionKey=injectionKey,
    )
    # Add a new injected_id column.
    output_catalog.add_column(col=list(range(len(output_catalog))), name="injected_id")
    # Remove unneccesary output columns.
    output_catalog.remove_column("patch")
    # Reorder columns
    mag_cols = [col for col in output_catalog.columns if f"_{col_mag}" in col]
    new_order = [
        "injected_id",
        "ra",
        "dec",
        "source_type",
        *mag_cols,
        "injection_id",
        "injection_draw_size",
        "injection_flag",
        "injected_isPatchInner",
        "injected_isTractInner",
        "injected_isPrimary",
    ]

    return output_catalog[new_order]


class ConsolidateInjectedCatalogsTask(PipelineTask):
    """Class for combining all tables in a collection of input catalogs
    into one table.
    """

    _DefaultName = "consolidateInjectedCatalogsTask"
    ConfigClass = ConsolidateInjectedCatalogsConfig

    def runQuantum(self, butlerQC, input_refs, output_refs):
        inputs = butlerQC.get(input_refs)
        catalog_dict, tract = _get_catalogs(inputs, input_refs)
        outputs = self.run(
            catalog_dict,
            inputs["skyMap"],
            tract,
            self.config.pixel_match_radius,
            self.config.get_catalogs_from_butler,
            self.config.col_ra,
            self.config.col_dec,
            self.config.col_mag,
            self.config.isPatchInnerKey,
            self.config.isTractInnerKey,
            self.config.isPrimaryKey,
            self.config.injectionKey,
        )
        butlerQC.put(outputs, output_refs)

    def run(
        self,
        catalog_dict: dict,
        skymap: BaseSkyMap,
        tract: int,
        pixel_match_radius: float = 0.1,
        get_catalogs_from_butler: bool = True,
        col_ra: str = "ra",
        col_dec: str = "dec",
        col_mag: str = "mag",
        isPatchInnerKey="injected_isPatchInner",
        isTractInnerKey="injected_isTractInner",
        isPrimaryKey="injected_isPrimary",
        injectionKey="injection_flag",
    ) -> Table:
        """Consolidate all tables in catalog_dict into one table.

        catalog_dict: `dict`
            A dictionary with photometric bands for keys and astropy tables for
            items.
        skymap: `lsst.skymap.BaseSkyMap`
            A base skymap.
        tract: `int`
            The tract where sources have been injected.
        pixel_match_radius: `float`
            Match radius in pixels to use for self-matching catalogs across
            different bands.
        col_ra: `str`
            Column name for right ascension (in degrees).
        col_dec: `str`
            Column name for declination (in degrees).
        col_mag: `str`
            Column name for magnitude.
        isPatchInnerKey: `str`
            Column name for the isPatchInner flag.
        isTractInnerKey: `str`
            Column name for the isTractInner flag.
        isPrimaryKey: `str`
            Column name for the isPrimary flag.
        injectionKey: `str`
            Column name for the injection flag.

        Returns
        -------
        output_struct : `lsst.pipe.base.Struct`
            contains :
                multiband_catalog: `astropy.table.Table`
                    A single table containing all information of the separate
                    tables in catalog_dict
        """
        output_catalog = consolidate_injected_deepCoadd_catalogs(
            catalog_dict,
            skymap,
            tract,
            pixel_match_radius,
            get_catalogs_from_butler,
            col_ra,
            col_dec,
            col_mag,
            isPatchInnerKey,
            isTractInnerKey,
            isPrimaryKey,
            injectionKey,
        )
        output_struct = Struct(output_catalog=output_catalog)
        return output_struct
