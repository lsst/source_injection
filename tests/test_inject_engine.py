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

import logging
import unittest
from types import GeneratorType

import galsim
import numpy as np
import pytest
from galsim import BoundsI, GSObject

import lsst.utils.tests
from lsst.geom import Point2D, SpherePoint, degrees
from lsst.images import YX
from lsst.images.cells import CellIJ
from lsst.source.injection.inject_engine import (
    generate_galsim_objects,
    get_gain_map,
    infer_gain_from_image,
    inject_galsim_objects_into_exposure,
    make_galsim_object,
)
from lsst.source.injection.utils.test_utils import (
    make_test_cell_coadd,
    make_test_exposure,
    make_test_injection_catalog,
)
from lsst.utils.tests import TestCase


class InjectEngineTestCase(TestCase):
    """Test the inject_engine.py module."""

    def setUp(self):
        """Set up synthetic source injection inputs.

        This method sets up a noisy synthetic image with Gaussian PSFs injected
        into the frame, an example source injection catalog, and a generator of
        GalSim objects intended for injection.
        """
        self.exposure = make_test_exposure()
        # Mark a cell as missing to test that it's ignored, although the data
        # in it may still actually be fine
        self.cell_bad = CellIJ(1, 1)
        self.cell_coadd = make_test_cell_coadd(
            exposure=self.exposure,
            cell_shape=YX(x=35, y=35),
            psf_shape=(33, 33),
            missing={self.cell_bad},
            band="r",
        )
        cen_cell_bad = self.exposure.wcs.pixelToSky(
            self.cell_coadd.grid.bbox_of(CellIJ(1, 1)).to_legacy().getCenter()
        )
        injection_catalog = make_test_injection_catalog(
            self.exposure.getWcs(),
            self.exposure.getBBox(),
        )
        row_last = injection_catalog[-1]
        injection_catalog.add_row(
            {
                "ra": cen_cell_bad.getRa().asDegrees(),
                "dec": cen_cell_bad.getDec().asDegrees(),
                "mag": row_last["mag"],
                "source_type": row_last["source_type"],
            }
        )
        # Add one galaxy that should span multiple cells
        self.cell_galaxy = CellIJ(1, 1)
        cen_cell_galaxy = self.exposure.wcs.pixelToSky(
            self.cell_coadd.grid.bbox_of(CellIJ(1, 1)).to_legacy().getCenter()
        )
        n_rows = len(injection_catalog)
        columns_sersic = ("n", "half_light_radius", "q", "beta")
        injection_catalog.add_columns(
            cols=tuple(
                np.ma.masked_array(data=np.zeros(n_rows, dtype=float), mask=np.ones(n_rows, dtype=bool))
                for _ in range(len(columns_sersic))
            ),
            names=columns_sersic,
        )
        injection_catalog.add_row(
            {
                "ra": cen_cell_galaxy.getRa().asDegrees(),
                "dec": cen_cell_galaxy.getDec().asDegrees(),
                "mag": row_last["mag"],
                "source_type": "Sersic",
                "n": 4.0,
                "half_light_radius": 20.0,
                "q": 0.8,
                "beta": 15.0,
            }
        )
        self.injection_catalog = injection_catalog

        self.galsim_objects = generate_galsim_objects(
            injection_catalog=self.injection_catalog,
            photo_calib=self.exposure.photoCalib,
            wcs=self.exposure.wcs,
            fits_alignment="wcs",
            stamp_prefix="",
        )
        self.photoCalib = self.exposure.getPhotoCalib()
        self.inst_fluxes = [
            float(self.photoCalib.magnitudeToInstFlux(mag)) for mag in self.injection_catalog["mag"]
        ]

    def tearDown(self):
        del self.exposure
        del self.cell_coadd
        del self.injection_catalog
        del self.galsim_objects
        del self.photoCalib
        del self.inst_fluxes

    def test_make_galsim_object(self):
        source_data = self.injection_catalog[0]
        sky_coords = SpherePoint(float(source_data["ra"]), float(source_data["dec"]), degrees)
        pixel_coords = self.exposure.wcs.skyToPixel(sky_coords)
        inst_flux = self.exposure.photoCalib.magnitudeToInstFlux(source_data["mag"], pixel_coords)
        object = make_galsim_object(
            source_data=source_data,
            source_type=source_data["source_type"],
            inst_flux=inst_flux,
        )
        self.assertIsInstance(object, GSObject)
        self.assertIsInstance(object, getattr(galsim, source_data["source_type"]))

    def test_generate_galsim_objects(self):
        self.assertTrue(isinstance(self.galsim_objects, GeneratorType))
        for galsim_object in self.galsim_objects:
            self.assertIsInstance(galsim_object, tuple)
            self.assertEqual(len(galsim_object), 4)
            self.assertIsInstance(galsim_object[0], SpherePoint)  # RA/Dec
            self.assertIsInstance(galsim_object[1], Point2D)  # x/y
            self.assertIsInstance(galsim_object[2], int)  # draw size
            self.assertIsInstance(galsim_object[3], GSObject)  # GSObject

    def test_infer_gain_nonexistent_mask(self):
        """Test that infer_gain_from_image returns a gain value and logs a
        warning when provided with nonexistent mask plane names.
        """
        logger = logging.getLogger(__name__)
        with self.assertLogs(logger, level="WARNING") as cm:
            gain = infer_gain_from_image(
                self.exposure,
                bad_mask_names=["NONEXISTENT_MASK1", "NONEXISTENT_MASK2"],
                logger=logger,
            )
        self.assertTrue(any("NONEXISTENT_MASK1" in msg for msg in cm.output))
        self.assertIsInstance(gain, float)
        self.assertTrue(np.isfinite(gain))
        gain_no_mask = infer_gain_from_image(self.exposure, bad_mask_names=[])
        self.assertAlmostEqual(gain, gain_no_mask)

    def test_infer_gain_no_valid_pixels(self):
        """Test that infer_gain_from_image returns NaN (rather than raising)
        when a region has no valid pixels to fit.
        """
        logger = logging.getLogger(__name__)
        # Every pixel flagged with a bad mask plane.
        all_bad = make_test_exposure()
        all_bad.mask.addMaskPlane("BAD")
        all_bad.mask.array[:] = all_bad.mask.getPlaneBitMask("BAD")
        with self.assertLogs(logger, level="WARNING"):
            gain = infer_gain_from_image(all_bad, bad_mask_names=["BAD"], logger=logger)
        self.assertTrue(np.isnan(gain))

        # Every variance pixel non-finite.
        all_nan = make_test_exposure()
        all_nan.variance.array[:] = np.nan
        self.assertTrue(np.isnan(infer_gain_from_image(all_nan, bad_mask_names=[])))

    def test_get_gain_map_no_valid_pixels(self):
        """Test that get_gain_map produces a finite, positive map even when no
        region can be fit, falling back to unit gain.
        """
        all_bad = make_test_exposure()
        all_bad.mask.addMaskPlane("BAD")
        all_bad.mask.array[:] = all_bad.mask.getPlaneBitMask("BAD")
        gain_map = get_gain_map(all_bad, bad_mask_names=["BAD"])
        self.assertTrue(np.all(np.isfinite(gain_map.array)))
        self.assertTrue(np.all(gain_map.array > 0))

    def test_inject_galsim_objects_into_cell_coadd(self):
        self._test_inject_galsim_objects_into_exposure(self.cell_coadd, True)

    def test_inject_galsim_objects_into_exposure(self):
        self._test_inject_galsim_objects_into_exposure(self.exposure, False)

    def _test_inject_galsim_objects_into_exposure(self, exposure, is_cell: bool = True):
        flux0 = np.sum(exposure.image.array)
        for injection_core_size in (None, 0, 1.5):
            with pytest.raises(ValueError):
                inject_galsim_objects_into_exposure(
                    exposure=exposure,
                    objects=(),
                    injection_core_size=injection_core_size,
                )
        injected_outputs = inject_galsim_objects_into_exposure(
            exposure=exposure,
            objects=self.galsim_objects,
            mask_plane_name="INJECTED",
            calib_flux_radius=12.0,
            draw_size_max=1000,
            add_noise=False,
            injection_core_size=5,
        )
        draw_sizes, common_bounds, fft_size_errors, psf_compute_errors = injected_outputs
        self.assertAlmostEqual(
            np.sum(exposure.image.array) - flux0,
            np.sum(self.inst_fluxes),
            delta=0.00015 * np.sum(self.inst_fluxes),
        )
        self.assertEqual(len(draw_sizes), len(self.injection_catalog["ra"]))
        self.assertTrue(all(isinstance(injected_output, list) for injected_output in injected_outputs))
        self.assertTrue(all(isinstance(item, int) for item in draw_sizes))
        self.assertTrue(all(isinstance(item, BoundsI) for item in common_bounds))  # common bounds
        self.assertTrue(all(isinstance(item, bool) for item in fft_size_errors))  # FFT size errors
        self.assertTrue(all(isinstance(item, bool) for item in psf_compute_errors))  # PSF compute errors
        mask_dict = exposure.mask.schema if is_cell else exposure.mask.getMaskPlaneDict()
        assert "INJECTED" in mask_dict
        assert "INJECTED_CORE" in mask_dict
        # Non-cell coadds shouldn't fail to inject here
        # Cell coadds actually should skip this object, but as it is,
        # the draw_size is the intended width/length of the injection box,
        # not the number of pixels that were actually injected into.
        assert draw_sizes[-1] > 0


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    """Test memory usage of functions in this script."""

    pass


def setup_module(module):
    """Configure pytest."""
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
