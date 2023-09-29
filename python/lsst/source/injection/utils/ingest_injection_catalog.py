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

__all__ = ["ingest_injection_catalog"]

import logging

import esutil
from astropy.table import Table, vstack
from lsst.daf.butler import Butler, CollectionType, DatasetRef, DatasetType, DimensionUniverse


def ingest_injection_catalog(
    writeable_butler: Butler,
    table: Table | list[Table],
    band: str,
    output_collection: str,
    dataset_type_name: str = "injection_catalog",
    log_level: int = logging.INFO,
) -> list[DatasetRef]:
    """Ingest a source table into the butler.

    This function ingests either a single astropy Table or list of Tables into
    the Butler. If a list of Tables is provided, these will be vertically
    stacked together into one single Table for ingestion. Input source tables
    are expected to contain the columns ``ra`` and ``dec``, with data in units
    of degrees. This spatial information will be used to shard up the source
    table on-disk using a depth 7 Hierarchical Triangular Mesh (HTM7) format.
    HTM7 trixels have an area of ~0.315 square degreees.

    Parameters
    ----------
    writeable_butler : `lsst.daf.butler.Butler`
        An instantiated writeable butler.
    table : `astropy.table.Table` or `list` [`astropy.table.Table`]
        Input source table(s). Requires columns ``ra`` and ``dec`` in degrees.
    band : `str`
        Band associated with the input source table(s).
    output_collection : `str`
        Name of the output collection to ingest the consolidated table into.
    dataset_type_name : `str`, optional
        Dataset type name for the ingested consolidated table.
    log_level : `int`, optional
        The log level to use for logging.

    Returns
    -------
    dataset_refs : `list` [`lsst.daf.butler.DatasetRef`]
        List containing the dataset refs for the ingested source table.
    """
    # Instantiate logger.
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # Concatenate inputs to a single Table if required.
    if isinstance(table, list):
        table = vstack(table)

    # Register the dataset type and collection.
    table_dataset_type = DatasetType(
        dataset_type_name,
        ("htm7", "band"),
        "ArrowAstropy",
        universe=DimensionUniverse(),
    )
    writeable_butler.registry.registerDatasetType(table_dataset_type)
    writeable_butler.registry.registerCollection(
        name=output_collection,
        type=CollectionType.RUN,
    )

    # Generate HTM trixel indices.
    htm = esutil.htm.HTM(7)
    htm_indices = htm.lookup_id(table["ra"], table["dec"])  # type: ignore

    # Loop over all unique trixel indices and ingest.
    dataset_refs: list[DatasetRef] = []
    for htm_index in set(htm_indices):
        table_subset = table[htm_indices == htm_index]
        dataset_ref = writeable_butler.put(
            table_subset,
            table_dataset_type,
            {"htm7": htm_index, "band": band},
            run=output_collection,
        )
        dataset_refs.append(dataset_ref)

    # Log results and return.
    logger.info(
        "Ingested %d %s band %s DatasetRef%s into the butler.",
        len(dataset_refs),
        band,
        dataset_type_name,
        "" if len(dataset_refs) == 1 else "s",
    )
    return dataset_refs
