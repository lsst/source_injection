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
from typing import cast

import numpy as np
from astropy.table import Table, vstack
from astropy.table.column import MaskedColumn
from lsst.geom import SpherePoint, degrees
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


class ConsolidateInjectedCatalogsTask(PipelineTask):
    """Class for combining all tables in a collection of input catalogs
    into one table.
    """

    _DefaultName = "consolidateInjectedCatalogsTask"
    ConfigClass = ConsolidateInjectedCatalogsConfig

    def runQuantum(self, butlerQC, input_refs, output_refs):
        inputs = butlerQC.get(input_refs)
        catalog_dict, tract = self._get_catalogs(inputs, input_refs)
        outputs = self.run(catalog_dict, inputs["skyMap"], tract)
        butlerQC.put(outputs, output_refs)

    def _get_catalogs(
        self,
        inputs: dict,
        input_refs: InputQuantizedConnection,
    ) -> tuple[dict, int]:
        """Organize input catalogs into a dictionary with photometry band
        keys.
        """
        catalog_dict: dict[str, list] = {}
        tracts = set()
        for ref, catalog in zip(input_refs.input_catalogs, inputs["input_catalogs"]):
            band = ref.dataId.band.name
            if band not in catalog_dict:
                catalog_dict[band] = []
            # load the patch number to check for patch overlap duplicates later
            catalog["patch"] = ref.dataId.patch.id
            catalog_dict[band].append(catalog)
            tracts.add(ref.dataId.tract.id)
        # vstack all the catalogs for each band
        for band, catalog_list in catalog_dict.items():
            catalog_dict[band] = vstack(catalog_list)
        if len(tracts) != 1:
            raise RuntimeError(f"len({tracts=}) != 1")
        return (catalog_dict, list(tracts)[0])

    def _remove_duplicates(
        self,
        catalog: Table,
        tractInfo,
    ) -> Table:
        """Remove patch overlap duplicates."""
        self.config = cast(ConsolidateInjectedCatalogsConfig, self.config)
        # TODO: add isPatchInner column
        duplicates = []
        for ind, row in enumerate(catalog):
            spherePoint = SpherePoint(row[self.config.col_ra] * degrees, row[self.config.col_dec] * degrees)
            # remove any sources not in the tract
            if tractInfo.contains(spherePoint):
                # get the patch number from the source's ra,dec
                patchInfo = tractInfo.findPatch(spherePoint)
                # check against the patch column
                if row["patch"] != patchInfo.sequential_index:
                    duplicates.append(ind)
            else:
                duplicates.append(ind)
        catalog.remove_rows(duplicates)
        return catalog

    def _make_multiband_catalog(
        self,
        bands: list,
        catalog_dict: dict,
        match_radius: float,
    ) -> Table:
        """Combine multiple band-specific catalogs into one multiband
        catalog.
        """
        self.config = cast(ConsolidateInjectedCatalogsConfig, self.config)
        # load the first catalog then loop to add info for the other bands
        multiband_catalog = catalog_dict[bands[0]]
        mag_suffix = self.config.col_mag
        multiband_catalog.rename_column(self.config.col_mag, f"{bands[0]}_{mag_suffix}")
        for band in bands[1:]:
            # make a column for the new band
            multiband_catalog.add_column([np.nan] * len(multiband_catalog), name=f"{band}_{mag_suffix}")
            # match the input catalog for this band to the existing
            # multiband catalog
            catalog_next_band = catalog_dict[band]
            catalog_next_band.rename_column(self.config.col_mag, f"{band}_{mag_suffix}")
            with Matcher(multiband_catalog[self.config.col_ra], multiband_catalog[self.config.col_dec]) as m:
                idx, multiband_match_inds, next_band_match_inds, dists = m.query_radius(
                    catalog_next_band[self.config.col_ra],
                    catalog_next_band[self.config.col_dec],
                    match_radius,
                    return_indices=True,
                )
            # if there are matches...
            if len(multiband_match_inds) > 0 and len(next_band_match_inds) > 0:
                # choose the coordinates in the brightest band
                for i, j in zip(multiband_match_inds, next_band_match_inds):
                    mags = []
                    for col in multiband_catalog.colnames:
                        if f"_{mag_suffix}" in col:
                            mags.append((col, multiband_catalog[i][col]))
                    bright_mag = min([x[1] for x in mags])
                    if catalog_next_band[f"{band}_{mag_suffix}"][j] < bright_mag:
                        multiband_catalog[self.config.col_ra][i] = catalog_next_band[self.config.col_ra][j]
                        multiband_catalog[self.config.col_dec][i] = catalog_next_band[self.config.col_dec][j]
                # TODO: Once multicomponent object support is added, make some
                # logic to pick the correct source_type.
                # Fill the new mag value.
                multiband_catalog[f"{band}_{mag_suffix}"][multiband_match_inds] = catalog_next_band[
                    f"{band}_{mag_suffix}"
                ][next_band_match_inds]
                # add rows for all the sources without matches
                not_next_band_match_inds = np.full(len(catalog_next_band), True, dtype=bool)
                not_next_band_match_inds[next_band_match_inds] = False
                multiband_catalog = vstack([multiband_catalog, catalog_next_band[not_next_band_match_inds]])
            # otherwise just stack the tables
            else:
                multiband_catalog = vstack([multiband_catalog, catalog_next_band])
        # fill any automatically masked values with NaNs
        if multiband_catalog.has_masked_columns:
            for col in multiband_catalog.columns:
                if isinstance(multiband_catalog[col], MaskedColumn):
                    multiband_catalog[col] = multiband_catalog[col].filled(np.nan)
        return multiband_catalog

    def run(
        self,
        catalog_dict: dict,
        skymap: BaseSkyMap,
        tract: int,
        pixel_match_radius: float = 0.1,
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

        Returns
        -------
        output_struct : `lsst.pipe.base.Struct`
            contains :
                multiband_catalog: `astropy.table.Table`
                    A single table containing all information of the separate
                    tables in catalog_dict
        """
        self.config = cast(ConsolidateInjectedCatalogsConfig, self.config)
        bands = list(catalog_dict.keys())
        # convert the pixel match radius to degrees
        tractInfo = skymap.generateTract(tract)
        tractWcs = tractInfo.getWcs()
        pixel_scale = tractWcs.getPixelScale()
        match_radius = pixel_match_radius * pixel_scale.asDegrees()
        for band in bands:
            catalog = catalog_dict[band]
            # remove patch overlap duplicates
            catalog = self._remove_duplicates(catalog, tractInfo)
        if len(bands) > 1:
            # match the catalogs across bands
            output_catalog = self._make_multiband_catalog(bands, catalog_dict, match_radius)
        else:
            output_catalog = catalog_dict[bands[0]]
            output_catalog.rename_column(self.config.col_mag, f"{bands[0]}_{self.config.col_mag}")
        # remove sources outside tract boundaries
        for index, (ra, dec) in enumerate(
            list(zip(output_catalog[self.config.col_ra], output_catalog[self.config.col_dec]))
        ):
            point = SpherePoint(ra * degrees, dec * degrees)
            if not tractInfo.contains(point):
                output_catalog.remove_row(index)
        # replace injection_id column with a new injected_id column
        output_catalog["injection_id"] = list(range(len(output_catalog)))
        output_catalog.rename_column("injection_id", "injected_id")
        output_struct = Struct(output_catalog=output_catalog)
        return output_struct
