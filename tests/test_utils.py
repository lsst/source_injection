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
import os
import unittest

import lsst.utils.tests
from lsst.daf.butler.tests import makeTestCollection, makeTestRepo
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.obs.base.instrument_tests import DummyCam
from lsst.pipe.base import Pipeline
from lsst.skymap.ringsSkyMap import RingsSkyMap, RingsSkyMapConfig
from lsst.source.injection import (
    ExposureInjectTask,
    consolidate_injected_deepCoadd_catalogs,
    ingest_injection_catalog,
    make_injection_pipeline,
)
from lsst.source.injection.utils.test_utils import (
    make_test_exposure,
    make_test_injection_catalog,
    make_test_reference_pipeline,
)
from lsst.utils.tests import MemoryTestCase, TestCase

TEST_DIR = os.path.abspath(os.path.dirname(__file__))


class SourceInjectionUtilsTestCase(TestCase):
    """Test the utility functions in the source_injection package."""

    @classmethod
    def setUpClass(cls):
        cls.root = makeTestTempDir(TEST_DIR)
        cls.creator_butler = makeTestRepo(cls.root)
        cls.writeable_butler = makeTestCollection(cls.creator_butler)
        # Register an instrument so we can get some bands.
        DummyCam().register(cls.writeable_butler.registry)
        skyMapConfig = RingsSkyMapConfig()
        skyMapConfig.numRings = 3
        cls.skyMap = RingsSkyMap(config=skyMapConfig)

    @classmethod
    def tearDownClass(cls):
        del cls.writeable_butler
        del cls.creator_butler
        del cls.skyMap
        removeTestTempDir(cls.root)

    def setUp(self):
        self.exposure = make_test_exposure()
        self.injection_catalog = make_test_injection_catalog(
            self.exposure.getWcs(),
            self.exposure.getBBox(),
        )
        self.reference_pipeline = make_test_reference_pipeline()
        self.injected_catalog = self.injection_catalog.copy()
        self.injected_catalog.add_columns(cols=[0, 0], names=["injection_draw_size", "injection_flag"])
        self.injected_catalog["injection_flag"][:5] = 1

    def tearDown(self):
        del self.exposure
        del self.injection_catalog
        del self.reference_pipeline
        del self.injected_catalog

    def test_generate_injection_catalog(self):
        self.assertEqual(len(self.injection_catalog), 30)
        expected_columns = {"injection_id", "ra", "dec", "source_type", "mag"}
        self.assertEqual(set(self.injection_catalog.columns), expected_columns)

    def test_make_injection_pipeline(self):
        injection_pipeline = Pipeline("reference_pipeline")
        injection_pipeline.addTask(ExposureInjectTask, "inject_exposure")
        merged_pipeline = make_injection_pipeline(
            dataset_type_name="postISRCCD",
            reference_pipeline=self.reference_pipeline,
            injection_pipeline=injection_pipeline,
            exclude_subsets=False,
            prefix="injected_",
            instrument="lsst.obs.subaru.HyperSuprimeCam",
            log_level=logging.DEBUG,
        )
        expanded_pipeline = merged_pipeline.toExpandedPipeline()
        expected_subset_tasks = ["isr", "inject_exposure", "characterizeImage"]
        merged_task_subsets = [merged_pipeline.findSubsetsWithLabel(x) for x in expected_subset_tasks]
        self.assertEqual(len(merged_task_subsets), len(expected_subset_tasks))
        for taskDef in expanded_pipeline:
            conns = taskDef.connections
            if taskDef.label == "isr":
                self.assertEqual(conns.outputExposure.name, "postISRCCD")
            elif taskDef.label == "inject_exposure":
                self.assertEqual(conns.input_exposure.name, "postISRCCD")
                self.assertEqual(conns.output_exposure.name, "injected_postISRCCD")
                self.assertEqual(conns.output_catalog.name, "injected_postISRCCD_catalog")
            elif taskDef.label == "characterizeImage":
                self.assertEqual(conns.exposure.name, "injected_postISRCCD")
                self.assertEqual(conns.characterized.name, "injected_icExp")
                self.assertEqual(conns.backgroundModel.name, "injected_icExpBackground")
                self.assertEqual(conns.sourceCat.name, "injected_icSrc")

    def test_ingest_injection_catalog(self):
        input_dataset_refs = ingest_injection_catalog(
            writeable_butler=self.writeable_butler,
            table=self.injection_catalog,
            band="g",
            output_collection="test_collection",
            dataset_type_name="injection_catalog",
            log_level=logging.DEBUG,
        )
        output_dataset_refs = self.writeable_butler.registry.queryDatasets(
            "injection_catalog",
            collections="test_collection",
        )
        self.assertEqual(len(input_dataset_refs), output_dataset_refs.count())
        input_ids = {x.id for x in input_dataset_refs}
        output_ids = {x.id for x in output_dataset_refs}
        self.assertEqual(input_ids, output_ids)
        injected_catalog = self.writeable_butler.get(input_dataset_refs[0])
        self.assertTrue(all(self.injection_catalog == injected_catalog))

    def test_consolidate_injected_catalogs(self):
        catalog_dict = {"g": self.injected_catalog, "r": self.injected_catalog}
        output_catalog = consolidate_injected_deepCoadd_catalogs(
            catalog_dict=catalog_dict,
            skymap=self.skyMap,
            tract=9,
            pixel_match_radius=0.1,
            get_catalogs_from_butler=False,
            col_ra="ra",
            col_dec="dec",
            col_mag="mag",
            isPatchInnerKey="injected_isPatchInner",
            isTractInnerKey="injected_isTractInner",
            isPrimaryKey="injected_isPrimary",
            injectionKey="injection_flag",
        )
        self.assertEqual(len(output_catalog), 30)
        expected_columns = {
            "injection_id",
            "injected_id",
            "ra",
            "dec",
            "source_type",
            "g_mag",
            "r_mag",
            "injection_draw_size",
            "injection_flag",
            "injected_isPatchInner",
            "injected_isTractInner",
            "injected_isPrimary",
        }
        self.assertEqual(set(output_catalog.columns), expected_columns)
        self.assertEqual(sum(output_catalog["injection_flag"]), 5)
        self.assertEqual(sum(output_catalog["injected_isPatchInner"]), 30)
        self.assertEqual(sum(output_catalog["injected_isTractInner"]), 30)
        self.assertEqual(sum(output_catalog["injected_isPrimary"]), 25)


class MemoryTestCase(MemoryTestCase):
    """Test memory usage of functions in this script."""

    pass


def setup_module(module):
    """Configure pytest."""
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
