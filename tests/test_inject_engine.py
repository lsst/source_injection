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
from galsim import BoundsI, GSObject

import lsst.utils.tests
from lsst.geom import Point2D, SpherePoint, degrees
from lsst.source.injection.inject_engine import (
    generate_galsim_objects,
    get_gain_map,
    infer_gain_from_image,
    inject_galsim_objects_into_exposure,
    make_galsim_object,
)
from lsst.source.injection.utils.test_utils import make_test_exposure, make_test_injection_catalog
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
        self.injection_catalog = make_test_injection_catalog(
            self.exposure.getWcs(),
            self.exposure.getBBox(),
        )
        self.galsim_objects = generate_galsim_objects(
            injection_catalog=self.injection_catalog,
            photo_calib=self.exposure.photoCalib,
            wcs=self.exposure.wcs,
            fits_alignment="wcs",
            stamp_prefix="",
        )

    def tearDown(self):
        del self.exposure
        del self.injection_catalog
        del self.galsim_objects

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

    def test_inject_galsim_objects_into_exposure(self):
        flux0 = np.sum(self.exposure.image.array)
        injected_outputs = inject_galsim_objects_into_exposure(
            exposure=self.exposure,
            objects=self.galsim_objects,
            mask_plane_name="INJECTED",
            calib_flux_radius=12.0,
            draw_size_max=1000,
            add_noise=False,
        )
        pc = self.exposure.getPhotoCalib()
        inst_fluxes = [float(pc.magnitudeToInstFlux(mag)) for mag in self.injection_catalog["mag"]]
        self.assertAlmostEqual(
            np.sum(self.exposure.image.array) - flux0,
            np.sum(inst_fluxes),
            delta=0.00015 * np.sum(inst_fluxes),
        )
        self.assertEqual(len(injected_outputs[0]), len(self.injection_catalog["ra"]))
        self.assertTrue(all(isinstance(injected_output, list) for injected_output in injected_outputs))
        self.assertTrue(all(isinstance(item, int) for item in injected_outputs[0]))  # draw sizes
        self.assertTrue(all(isinstance(item, BoundsI) for item in injected_outputs[1]))  # common bounds
        self.assertTrue(all(isinstance(item, bool) for item in injected_outputs[2]))  # FFT size errors
        self.assertTrue(all(isinstance(item, bool) for item in injected_outputs[3]))  # PSF compute errors


class MemoryTestCase(lsst.utils.tests.MemoryTestCase):
    """Test memory usage of functions in this script."""

    pass


def setup_module(module):
    """Configure pytest."""
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
