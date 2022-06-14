# This file is part of source injection
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Inject sources into visit-level imaging.
"""

import lsst.pipe.base.connectionTypes as cT
from . import inject_base

__all__ = ["VisitInjectConfig", "VisitInjectTask"]


class VisitInjectConnections(inject_base.BaseInjectConnections,
                             dimensions=("instrument", "visit", "detector"),
                             ):

    inputExposure = cT.Input(
        doc="Exposure to inject synthetic sources into.",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector")
    )

    outputExposure = cT.Output(
        doc="Injected Exposure.",
        name="{siPrefix}calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector")
    )


class VisitInjectConfig(inject_base.BaseInjectConfig,
                        pipelineConnections=VisitInjectConnections,
                        ):

    pass


class VisitInjectTask(inject_base.BaseInjectTask):
    """Class to inject into visit level images.
    """

    _DefaultName = "visitInjectTask"
    ConfigClass = VisitInjectConfig

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        detectorId = inputs["inputExposure"].getDetector().getId()

        expWcs = inputs["inputExposure"].getWcs()
        tractId = inputs["skyMap"].findTract(
            expWcs.pixelToSky(inputs["inputExposure"].getBBox().getCenter())).tract_id
        if not self.config.doApplyExternalGlobalSkyWcs and not self.config.doApplyExternalTractSkyWcs:
            inputs["wcs"] = expWcs
        elif self.config.doApplyExternalGlobalSkyWcs:
            externalSkyWcsCatalog = inputs["externalSkyWcsGlobalCatalog"]
            row = externalSkyWcsCatalog.find(detectorId)
            inputs["wcs"] = row.getWcs()
        elif self.config.doApplyExternalTractSkyWcs:
            externalSkyWcsCatalogList = inputs["externalSkyWcsTractCatalog"]
            externalSkyWcsCatalog = None
            for externalSkyWcsCatalogRef in externalSkyWcsCatalogList:
                if externalSkyWcsCatalogRef.dataId["tract"] == tractId:
                    externalSkyWcsCatalog = externalSkyWcsCatalogRef.get(
                        datasetType=self.config.connections.externalSkyWcsTractCatalog)
                    break
            if externalSkyWcsCatalog is None:
                usedTract = externalSkyWcsCatalogList[-1].dataId["tract"]
                self.log.warn(
                    f"Warning, external SkyWcs for tract {tractId} not found. Using tract {usedTract} "
                    "instead.")
                externalSkyWcsCatalog = externalSkyWcsCatalogList[-1].get(
                    datasetType=self.config.connections.externalSkyWcsTractCatalog)
            row = externalSkyWcsCatalog.find(detectorId)
            inputs["wcs"] = row.getWcs()

        if not self.config.doApplyExternalGlobalPhotoCalib and not self.config.doApplyExternalTractPhotoCalib:
            inputs["photoCalib"] = inputs["inputExposure"].getPhotoCalib()
        elif self.config.doApplyExternalGlobalPhotoCalib:
            externalPhotoCalibCatalog = inputs["externalPhotoCalibGlobalCatalog"]
            row = externalPhotoCalibCatalog.find(detectorId)
            inputs["photoCalib"] = row.getPhotoCalib()
        elif self.config.doApplyExternalTractPhotoCalib:
            externalPhotoCalibCatalogList = inputs["externalPhotoCalibTractCatalog"]
            externalPhotoCalibCatalog = None
            for externalPhotoCalibCatalogRef in externalPhotoCalibCatalogList:
                if externalPhotoCalibCatalogRef.dataId["tract"] == tractId:
                    externalPhotoCalibCatalog = externalPhotoCalibCatalogRef.get(
                        datasetType=self.config.connections.externalPhotoCalibTractCatalog)
                    break
            if externalPhotoCalibCatalog is None:
                usedTract = externalPhotoCalibCatalogList[-1].dataId["tract"]
                self.log.warn(
                    f"Warning, external PhotoCalib for tract {tractId} not found. Using tract {usedTract} "
                    "instead.")
                externalPhotoCalibCatalog = externalPhotoCalibCatalogList[-1].get(
                    datasetType=self.config.connections.externalPhotoCalibTractCatalog)
            row = externalPhotoCalibCatalog.find(detectorId)
            inputs["photoCalib"] = row.getPhotoCalib()

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)
