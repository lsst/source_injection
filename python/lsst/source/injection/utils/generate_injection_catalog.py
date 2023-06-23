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

__all__ = ["generate_injection_catalog"]

import hashlib
import itertools
import logging
from collections.abc import Sequence
from typing import Any

import numpy as np
from astropy.table import Table, hstack
from lsst.afw.geom import SkyWcs
from scipy.stats import qmc


def generate_injection_catalog(
    ra_lim: Sequence[float],
    dec_lim: Sequence[float],
    wcs: SkyWcs = None,
    number: int = 1,
    density: int | None = None,
    seed: Any = None,
    log_level: int = logging.INFO,
    **kwargs: Any,
) -> Table:
    """Generate a synthetic source injection catalog.

    This function generates a synthetic source injection catalog from user
    supplied input parameters. The catalog is returned as an astropy Table.

    On-sky source positions are generated using the quasi-random Halton
    sequence. By default, the Halton sequence is seeded using the product of
    the right ascension and declination limit ranges. This ensures that the
    same sequence is always generated for the same limits. This seed may be
    overridden by specifying the ``seed`` parameter.

    A unique injection ID is generated for each source. The injection ID
    encodes two pieces of information: the unique source identification number
    and the version number of the source as specified by the ``number``
    parameter. To achieve this, the unique source ID number is multiplied by
    `10**n` such that the sum of the multiplied source ID number with the
    unique repeated version number will always be unique. For example, an
    injection catalog with `number = 3` versions of each source will have
    injection IDs: 0, 1, 2, 10, 11, 12, 20, 21, 22, etc. If `number = 20`, then
    the injection IDs will be: 0, 1, 2, ..., 17, 18, 19, 100, 101, 102, etc.
    If `number = 1` (default) then the injection ID will be a simple sequential
    list of integers.

    Parameters
    ----------
    ra_lim : `Sequence` [`float`]
        The right ascension limits of the catalog in degrees.
    dec_lim : `Sequence` [`float`]
        The declination limits of the catalog in degrees.
    wcs : `lsst.afw.geom.SkyWcs`, optional
        The WCS associated with these data. If not given or ``None`` (default),
        the catalog will be generated using Cartesian geometry.
    number : `int`, optional
        The number of times to generate each unique combination of input
        parameters. The default is 1 (i.e., no repeats). This will be ignored
        if ``density`` is specified.
    density : `int` | None, optional
        The desired source density in sources per square degree. If given, the
        ``number`` parameter will be ignored. Instead, the number of unique
        parameter combination generations will be calculated to achieve the
        desired density. The default is `None` (i.e., no density calculation).
    seed : `Any`, optional
        The seed to use for the Halton sequence. If not given or ``None``
        (default), the seed will be set using the product of the right
        ascension and declination limit ranges.
    log_level : `int`, optional
        The log level to use for logging.
    **kwargs : `Any`
        The input parameters used to generate the catalog. Each parameter key
        will be used as a column name in the catalog. The values are the unique
        values for that parameter. The output catalog will contain a row for
        each unique combination of input parameters and be generated the number
        of times specified by ``number``.

    Returns
    -------
    table : `astropy.table.Table`
        The fully populated synthetic source injection catalog. The catalog
        will contain an automatically generated ``injection_id`` column that
        is unique for each source. The injection ID encodes two pieces of
        information: the unique source identification number and the repeated
        version number of the source as defined by the ``number`` parameter.
    """
    # Instantiate logger.
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # Parse optional keyword input parameters.
    values: list[Any] = [np.atleast_1d(x) for x in kwargs.values()]

    # Determine the BBox limits and pixel scale.
    if wcs:
        sky_corners = list(itertools.product(ra_lim, dec_lim))
        ra_corners, dec_corners = np.array(sky_corners).T
        x_corners, y_corners = wcs.skyToPixelArray(ra_corners, dec_corners, degrees=True)
        xlim: Any = np.percentile(x_corners, [0, 100])
        ylim: Any = np.percentile(y_corners, [0, 100])
    else:
        xlim, ylim = ra_lim, dec_lim

    # Automatically calculate the number of generations if density is given.
    if density:
        dec_lim_rad = np.deg2rad(dec_lim)
        area = ((180 / np.pi) * np.diff(ra_lim) * (np.sin(dec_lim_rad[1]) - np.sin(dec_lim_rad[0])))[0]
        rows = list(itertools.product(*values))
        native_density = len(rows) / area
        number = np.round(density / native_density).astype(int)
        if number > 0:
            logger.info(
                "Setting number of generations to %s, equivalent to %.1f sources per square degree.",
                number,
                number * native_density,
            )
        else:
            logger.warning("Requested source density would require number < 1; setting number = 1.")
            number = 1

    # Generate the fully expanded parameter table.
    values.append(range(number))
    keys = list(kwargs.keys())
    keys.append("version_id")
    param_table = Table(rows=list(itertools.product(*values)), names=keys)

    # Generate on-sky coordinate pairs.
    if not seed:
        seed = str(np.diff(ra_lim)[0] * np.diff(dec_lim)[0])
    # Random seed is the lower 32 bits of the hashed name.
    # We use hashlib.sha256 for guaranteed repeatability.
    hex_hash = hashlib.sha256(seed.encode("UTF-8")).hexdigest()
    hashed_seed = int("0x" + hex_hash, 0) & 0xFFFFFFFF
    sampler = qmc.Halton(d=2, seed=hashed_seed)
    sample = sampler.random(n=len(param_table))
    # Flip RA values if no WCS given.
    if not wcs:
        sample[:, 0] = 1 - sample[:, 0]
    xy_coords = Table(qmc.scale(sample, [xlim[0], ylim[0]], [xlim[1], ylim[1]]), names=("x", "y"))
    if wcs:
        ra_coords, dec_coords = wcs.pixelToSkyArray(xy_coords["x"], xy_coords["y"], degrees=True)
        sky_coords = Table([ra_coords, dec_coords], names=("ra", "dec"))
    else:
        sky_coords = Table(xy_coords, names=("ra", "dec"))
    # Perform an additional random permutation of the sky coordinate pairs to
    # minimize the potential for on-sky parameter correlations.
    rng = np.random.default_rng(hashed_seed)
    sky_coords = Table(rng.permutation(sky_coords))

    # Generate the unique injection ID and construct the final table.
    source_id = np.concatenate([([i] * number) for i in range(int(len(param_table) / number))])
    injection_id = param_table["version_id"] + source_id * int(10 ** np.ceil(np.log10(number)))
    injection_id.name = "injection_id"
    table = hstack([injection_id, sky_coords, param_table])
    table.remove_column("version_id")

    # Final logger report and return.
    if number == 1:
        extra_info = f"{len(table)} unique sources."
    else:
        extra_info = f"{len(table)} sources: {int(len(table)/number)} combinations repeated {number} times."
    logger.info("Generated an injection catalog containing %s", extra_info)
    return table
