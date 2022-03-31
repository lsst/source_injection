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
Inject sources into imaging.
"""
import galsim
import numpy as np
import pandas as pd

from lsst.pipe.base import PipelineTask, PipelineTaskConfig, PipelineTaskConnections

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import lsst.pipe.base.connectionTypes as cT
from lsst.skymap import BaseSkyMap
from lsst.geom import Box2D
from . import inject_engine

__all__ = ["BaseInjectConnections", "BaseInjectConfig", "BaseInjectTask"]


class BaseInjectConnections(PipelineTaskConnections,
                            dimensions=("instrument", "visit", "detector"),
                            defaultTemplates={"photoCalibName": "jointcal",
                                              "wcsName": "jointcal",
                                              "siPrefix": "si_",
                                              },
                            ):

    skyMap = cT.Input(
        doc="Input definition of geometry/bbox and projection/wcs for "
        "template exposures. Needed to test which tract to generate.",
        name=BaseSkyMap.SKYMAP_DATASET_TYPE_NAME,
        dimensions=("skymap",),
        storageClass="SkyMap",
    )

    injectionCatalogs = cT.Input(
        doc="Set of catalogs of sources to draw inputs from. We "
            "concatenate the tract catalogs for detectorVisits that cover "
            "multiple tracts.",
        name="{siPrefix}cat",
        storageClass="DataFrame",
        dimensions=("tract", "skymap"),
        deferLoad=True,
        multiple=True,
    )

    externalPhotoCalibTractCatalog = cT.Input(
        doc=("Per-tract, per-visit photometric calibrations. These catalogs use the "
             "detector id for the catalog id, sorted on id for fast lookup."),
        name="{photoCalibName}PhotoCalibCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit", "tract"),
        deferLoad=True,
        multiple=True,
    )

    externalPhotoCalibGlobalCatalog = cT.Input(
        doc=("Per-visit photometric calibrations. These catalogs use the "
             "detector id for the catalog id, sorted on id for fast lookup."),
        name="{photoCalibName}PhotoCalibCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
    )

    externalSkyWcsTractCatalog = cT.Input(
        doc=("Per-tract, per-visit wcs calibrations. These catalogs use the detector "
             "id for the catalog id, sorted on id for fast lookup."),
        name="{wcsName}SkyWcsCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit", "tract", "skymap"),
        deferLoad=True,
        multiple=True,
    )

    externalSkyWcsGlobalCatalog = cT.Input(
        doc=("Per-visit wcs calibrations computed globally (with no tract information). "
             "These catalogs use the detector id for the catalog id, sorted on id for "
             "fast lookup."),
        name="{wcsName}SkyWcsCatalog",
        storageClass="ExposureCatalog",
        dimensions=("instrument", "visit"),
    )

    def __init__(self, *, config=None):
        super().__init__(config=config)

        if not config.doApplyExternalTractPhotoCalib:
            self.inputs.remove("externalPhotoCalibTractCatalog")
        if not config.doApplyExternalGlobalPhotoCalib:
            self.inputs.remove("externalPhotoCalibGlobalCatalog")

        if not config.doApplyExternalTractSkyWcs:
            self.inputs.remove("externalSkyWcsTractCatalog")
        if not config.doApplyExternalGlobalSkyWcs:
            self.inputs.remove("externalSkyWcsGlobalCatalog")


class BaseInjectConfig(PipelineTaskConfig,
                       pipelineConnections=BaseInjectConnections,
                       ):

    # PhotoCalib

    doApplyExternalTractPhotoCalib = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Whether to apply an external photometric calibration via an "
            "`lsst.afw.image.PhotoCalib` object. Uses the "
            "`externalPhotoCalibName` config option to determine which "
            "calibration to use. Uses a per tract calibration."
    )

    doApplyExternalGlobalPhotoCalib = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Whether to apply an external photometric calibration via an "
            "`lsst.afw.image.PhotoCalib` object. Uses the "
            "`externalPhotoCalibName` config option to determine which "
            "calibration to use. Uses a global calibration."
    )

    externalPhotoCalibName = pexConfig.ChoiceField(
        doc="What type of external photo calib to use.",
        dtype=str,
        default="jointcal",
        allowed={"jointcal": "Use jointcal_photoCalib",
                 "fgcm": "Use fgcm_photoCalib",
                 "fgcm_tract": "Use fgcm_tract_photoCalib"}
    )

    # WCS

    doApplyExternalTractSkyWcs = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Whether to apply an external astrometric calibration via an "
            "`lsst.afw.geom.SkyWcs` object. Uses the "
            "`externalSkyWcsName` config option to determine which "
            "calibration to use. Uses a per tract calibration."
    )

    doApplyExternalGlobalSkyWcs = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Whether to apply an external astrometric calibration via an "
            "`lsst.afw.geom.SkyWcs` object. Uses the "
            "`externalSkyWcsName` config option to determine which "
            "calibration to use. Uses a global calibration."
    )

    externalSkyWcsName = pexConfig.ChoiceField(
        doc="What type of updated WCS calib to use.",
        dtype=str,
        default="jointcal",
        allowed={"jointcal": "Use jointcal_wcs"}
    )

    # source definitions and catalog manipulation

    doMatchVisit = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="Match visit to trim the input SSI catalog"
    )

    doCleanCat = pexConfig.Field(
        doc="If true removes bad sources from the catalog.",
        dtype=bool,
        default=True,
    )

    calibFluxRadius = pexConfig.Field(
        doc="Aperture radius (in pixels) that was used to define the calibration for this image+catalog. "
        "This will be used to produce the correct instrumental fluxes within the radius. "
        "This value should match that of the field defined in slot_CalibFlux_instFlux.",
        dtype=float,
        default=9.0,  # NOTE: this should be 12.0, awaiting and RC2 run with DM-34698
    )

    doSubSelectSources = pexConfig.Field(
        doc="Set to True if you wish to sub select sources to be input based on the value in the column"
            "set in the sourceSelectionColName config option.",
        dtype=bool,
        default=False
    )

    injectImages = pexConfig.Field(
        doc="Inject images directly? True or False.",
        dtype=bool,
        default=False,
    )

    doProcessAllDataIds = pexConfig.Field(
        doc="If True, all input data IDs will be processed, even those containing no synthetic sources.",
        dtype=bool,
        default=False,
    )

    trimBuffer = pexConfig.Field(
        doc="Size of the pixel buffer surrounding the image. Only those syntetic sources with a centroid"
        "falling within the image+buffer region will be considered for source injection.",
        dtype=int,
        default=100,
    )

    sourceType = pexConfig.Field(
        doc="The column name for the source type used in the synthetic source catalog.",
        dtype=str,
        default="sourceType",
    )

    fits_alignment = pexConfig.ChoiceField(
        doc="How should injections from FITS files be aligned?",
        dtype=str,
        allowed={
            "wcs": (
                "Input image will be transformed such that the local WCS in "
                "the FITS header matches the local WCS in the target image. "
                "I.e., North, East, and angular distances in the input image "
                "will match North, East, and angular distances in the target "
                "image."
            ),
            "pixel": (
                "Input image will _not_ be transformed.  Up, right, and pixel "
                "distances in the input image will match up, right and pixel "
                "distances in the target image."
            )
        },
        default="pixel"
    )

    # Input catalog config variables

    ra_col = pexConfig.Field(
        doc="Source catalog column name for RA (in radians).",
        dtype=str,
        default="ra",
    )

    dec_col = pexConfig.Field(
        doc="Source catalog column name for dec (in radians).",
        dtype=str,
        default="dec",
    )

    mag_col = pexConfig.Field(
        doc="Source catalog column name template for magnitudes, in the format "
            "``filter name``_mag_col.  E.g., if this config variable is set to "
            "``%s_mag``, then the i-band magnitude will be searched for in the "
            "``i_mag`` column of the source catalog.",
        dtype=str,
        default="%s_mag"
    )

    bulge_semimajor_col = pexConfig.Field(
        doc="Source catalog column name for the semimajor axis (in arcseconds) "
            "of the bulge half-light ellipse.",
        dtype=str,
        default="bulge_semimajor",
    )

    bulge_axis_ratio_col = pexConfig.Field(
        doc="Source catalog column name for the axis ratio of the bulge "
            "half-light ellipse.",
        dtype=str,
        default="bulge_axis_ratio",
    )

    bulge_pa_col = pexConfig.Field(
        doc="Source catalog column name for the position angle (measured from "
            "North through East in degrees) of the semimajor axis of the bulge "
            "half-light ellipse.",
        dtype=str,
        default="bulge_pa",
    )

    bulge_n_col = pexConfig.Field(
        doc="Source catalog column name for the Sersic index of the bulge.",
        dtype=str,
        default="bulge_n",
    )

    disk_semimajor_col = pexConfig.Field(
        doc="Source catalog column name for the semimajor axis (in arcseconds) "
            "of the disk half-light ellipse.",
        dtype=str,
        default="disk_semimajor",
    )

    disk_axis_ratio_col = pexConfig.Field(
        doc="Source catalog column name for the axis ratio of the disk "
            "half-light ellipse.",
        dtype=str,
        default="disk_axis_ratio",
    )

    disk_pa_col = pexConfig.Field(
        doc="Source catalog column name for the position angle (measured from "
            "North through East in degrees) of the semimajor axis of the disk "
            "half-light ellipse.",
        dtype=str,
        default="disk_pa",
    )

    disk_n_col = pexConfig.Field(
        doc="Source catalog column name for the Sersic index of the disk.",
        dtype=str,
        default="disk_n",
    )

    bulge_disk_flux_ratio_col = pexConfig.Field(
        doc="Source catalog column name for the bulge/disk flux ratio.",
        dtype=str,
        default="bulge_disk_flux_ratio",
    )

    select_col = pexConfig.Field(
        doc="Source catalog column name to be used to select which sources to "
            "add.",
        dtype=str,
        default="select",
    )

    def setDefaults(self):
        super().setDefaults()


class BaseInjectTask(PipelineTask):
    """Class to inject into visit level images.
    """

    _DefaultName = "baseInjectTask"
    ConfigClass = BaseInjectConfig

    def runQuantum(self, butlerQC, inputRefs, outputRefs):
        inputs = butlerQC.get(inputRefs)
        inputs["wcs"] = inputs["inputExposure"].getWcs()
        inputs["photoCalib"] = inputs["inputExposure"].getPhotoCalib()

        outputs = self.run(**inputs)
        butlerQC.put(outputs, outputRefs)

    def run(self, injectionCatalogs, inputExposure, skyMap, wcs, photoCalib, **kwargs):
        """Inject sources into a calexp and then run detection, deblending and
        measurement.

        Parameters
        ----------
        injectionCatalogs : `list` of `lsst.daf.butler.DeferredDatasetHandle`
            Set of tract level synthetic catalogs that potentially cover this
            detectorVisit.
        inputExposure : `lsst.afw.image.exposure.ExposureF`
            The exposure to inject the synthetic sources into.
        skyMap : `lsst.skymap.SkyMap`
            SkyMap defining the tracts and patches the synthetic sournces are
            stored over.
        wcs : `lsst.afw.geom.SkyWcs`
            WCS to use to inject synthetic sources.
        photoCalib : `lsst.afw.image.photoCalib.PhotoCalib`
            Photometric calibration to be used to calibrate the synthetic
            sources.
        exposureIdInfo : `lsst.obs.base.ExposureIdInfo`
        sfdSourceCat : `lsst.afw.table.SourceCatalog`
            Default : None
            Catalog produced by singleFrameDriver, needed to copy some
            calibration flags from.

        Returns
        -------
        resultStruct : `lsst.pipe.base.struct.Struct`
            contains : outputExposure : `lsst.afw.image.exposure.ExposureF`
                       outputCat : `lsst.afw.table.source.SourceCatalog`

        Notes
        -----
        Adds pixel coordinates for each source to the siCat and removes
        sources with bulge or disk half light radius = 0
        (if ``config.cleanCat = True``). These columns are called ``x`` and
        ``y`` and are in pixels.

        Adds the ``SI`` mask plane to the exposure which is then set by
        `inject_sources` to mark where synthetic sources have been injected.
        Uses the information in the ``siCat`` to make synthetic galaxies
        using GalSim and stars using the PSF models from the PSF information
        for the calexp. These are then added to the calexp, returning the
        calexp with syntetics included.

        GalSim galaxies are made using a double sersic profile, one for the
        bulge and one for the disk. This object is then convolved with the PSF.

        If exposureIdInfo is not provided then SourceCatalog IDs will not be
        globally unique.
        """
        siCat = self.composeSiCat(injectionCatalogs, skyMap)

        if self.config.doMatchVisit:
            # Trim the siCat to select a particular visit
            siCat = siCat[inputExposure.getInfo().getVisitInfo().getId() == siCat["visit"]]

        # Attach overriding wcs and photoCalib to inputExposure; keep originals
        # so we can reset at the end.
        origWcs = inputExposure.getWcs()
        origPhotoCalib = inputExposure.getPhotoCalib()
        inputExposure.setWcs(wcs)
        inputExposure.setPhotoCalib(photoCalib)

        band = inputExposure.getFilterLabel().bandLabel
        siCat = self._standardizeColumns(siCat, band)

        siCat = inject_engine.addPixCoords(siCat, inputExposure)
        siCat = self.trimSiCat(siCat, inputExposure)

        if len(siCat) > 0:
            if not self.config.injectImages:
                if isinstance(siCat[self.config.sourceType].iloc[0], str):
                    galCheckVal = "galaxy"
                    starCheckVal = "star"
                elif isinstance(siCat[self.config.sourceType].iloc[0], bytes):
                    galCheckVal = b"galaxy"
                    starCheckVal = b"star"
                elif isinstance(siCat[self.config.sourceType].iloc[0], (int, float)):
                    galCheckVal = 1
                    starCheckVal = 0
                else:
                    raise TypeError("sourceType column does not have required type: should be str, bytes or "
                                    "int/float.")
                if self.config.doCleanCat:
                    siCat = self.cleanCat(siCat, starCheckVal)

                generator = inject_engine._generateGSObjectsFromCatalog(
                    inputExposure, siCat, galCheckVal, starCheckVal, sourceTypeColumn=self.config.sourceType,
                    logger=self.log)
            else:
                generator = inject_engine._generateGSObjectsFromImages(inputExposure, siCat)
            inject_engine.inject_sources(inputExposure, generator,
                                         calibFluxRadius=self.config.calibFluxRadius, logger=self.log)
        elif len(siCat) == 0 and self.config.doProcessAllDataIds:
            self.log.warning("No synthetic sources found for this dataRef; processing anyway.")
            inputExposure.mask.addMaskPlane("SI")
        else:
            raise RuntimeError("No synthetic sources found for this dataRef.")

        # restore original inputExposure WCS and photoCalib
        inputExposure.setWcs(origWcs)
        inputExposure.setPhotoCalib(origPhotoCalib)

        resultStruct = pipeBase.Struct(outputExposure=inputExposure, outputCat=siCat)
        return resultStruct

    def composeSiCat(self, injectionCatalogs, skyMap):
        """Concatenate the injectionCatalogs from tracts that may cover the
        exposure.

        Parameters
        ----------
        injectionCatalogs : `list` of `lsst.daf.butler.DeferredDatasetHandle`
            Set of synthetic source catalogs to concatenate.
        skyMap : `lsst.skymap.SkyMap`
            SkyMap defining the geometry of the tracts and patches.

        Returns
        -------
        combinedSiCat : `pandas.DataFrame`
            All synthetic sources that cover the inner polygon of the tracts in
            this quantum.
        """
        if len(injectionCatalogs) == 1:
            return injectionCatalogs[0].get(
                datasetType=self.config.connections.injectionCatalogs)
        outputCat = []
        for siCatRef in injectionCatalogs:
            cat = siCatRef.get(
                datasetType=self.config.connections.injectionCatalogs)
            tractId = siCatRef.dataId["tract"]
            # Make sure all data is within the inner part of the tract.
            outputCat.append(cat[
                skyMap.findTractIdArray(cat[self.config.ra_col],
                                        cat[self.config.dec_col],
                                        degrees=False)
                == tractId])

        return pd.concat(outputCat)

    def _standardizeColumns(self, siCat, band):
        """Use config variables to 'standardize' the expected columns and
        column names in the input catalog. This method replaces all column
        names in the config with hard-coded internal names.

        Parameters
        ----------
        siCat : `pandas.core.frame.DataFrame`
            The catalog of synthetic sources to be injected.
        band : `str`
            Label for the current band being processed.

        Returns
        -------
        outCat : `pandas.core.frame.DataFrame`
            The standardized catalog of synthetic sources.
        """
        cfg = self.config
        replace_dict = {}

        def add_to_replace_dict(new_name, std_name):
            if new_name in siCat.columns:
                replace_dict[new_name] = std_name
            else:
                raise ValueError(f"Could not determine column for {std_name}.")

        # RA, Dec and mag are always required. Do these first.
        for new_name, std_name in [(cfg.ra_col, "ra"),
                                   (cfg.dec_col, "dec"),
                                   (cfg.mag_col%band, "mag"),
                                   ]:
            add_to_replace_dict(new_name, std_name)

        # add default bulge-to-disk flux ratio column if none already specified
        if "bulge_disk_flux_ratio" not in siCat.columns:
            siCat["bulge_disk_flux_ratio"] = 1.0

        # Only handle GalSim params if not injecting images
        if not cfg.injectImages:
            for new_name, std_name in [(cfg.bulge_semimajor_col, "bulge_semimajor"),
                                       (cfg.bulge_axis_ratio_col, "bulge_axis_ratio"),
                                       (cfg.bulge_pa_col, "bulge_pa"),
                                       (cfg.bulge_n_col, "bulge_n"),
                                       (cfg.disk_semimajor_col, "disk_semimajor"),
                                       (cfg.disk_axis_ratio_col, "disk_axis_ratio"),
                                       (cfg.disk_pa_col, "disk_pa"),
                                       (cfg.disk_n_col, "disk_n"),
                                       (cfg.bulge_disk_flux_ratio_col, "bulge_disk_flux_ratio"),
                                       ]:
                add_to_replace_dict(new_name, std_name)

        # Only handle flag column selection if specifically requested
        if cfg.doSubSelectSources:
            add_to_replace_dict(cfg.select_col, "select")

        # standardize columns
        siCat = siCat.rename(columns=replace_dict, copy=False)

        return siCat

    def trimSiCat(self, siCat, image):
        """Trim the synthetic catalog to about the size of the input image.

        `siCat` must be processed with addPixCoords before using this method.

        Parameters
        ----------
        siCat : `pandas.core.frame.DataFrame`
            The catalog of synthetic sources to be injected.
        image : `lsst.afw.image.exposure.ExposureF`
            The image into which the synthetic sources should be injected.

        Returns
        -------
        siCat : `pandas.core.frame.DataFrame`
            The original siCat trimmed to the area of the image.
        """

        bbox = Box2D(image.getBBox()).dilatedBy(self.config.trimBuffer)
        xs = siCat["x"].values
        ys = siCat["y"].values

        isContained = xs >= bbox.minX
        isContained &= xs <= bbox.maxX
        isContained &= ys >= bbox.minY
        isContained &= ys <= bbox.maxY

        return siCat[isContained]

    def cleanCat(self, siCat, starCheckVal):
        """Remove rows from the synthetic catalog which have HLR = 0 for either
        the bulge or disk component. Also remove galaxies that have Sersic
        index outside the GalSim min and max allowed range (0.3 <= n <= 6.2).

        Parameters
        ----------
        siCat : `pandas.core.frame.DataFrame`
            The catalog of synthetic sources to be input
        starCheckVal : `str`, `bytes` or `int`
            The value set in the sourceType column to specifiy a star object.

        Returns
        -------
        siCat : `pandas.core.frame.DataFrame`
            The input catalog of synthetic sources with bad objects removed.
        """

        rowsToKeep = (((siCat["bulge_semimajor"] != 0.0) & (siCat["disk_semimajor"] != 0.0))
                      | (siCat[self.config.sourceType] == starCheckVal))
        numRowsNotUsed = len(siCat) - len(np.where(rowsToKeep)[0])
        self.log.info("Removing %d rows with HLR = 0 for either the bulge or disk", numRowsNotUsed)
        siCat = siCat[rowsToKeep]

        minN = galsim.Sersic._minimum_n
        maxN = galsim.Sersic._maximum_n
        rowsWithGoodSersic = (((siCat["bulge_n"] >= minN) & (siCat["bulge_n"] <= maxN)
                              & (siCat["disk_n"] >= minN) & (siCat["disk_n"] <= maxN))
                              | (siCat[self.config.sourceType] == starCheckVal))
        numRowsNotUsed = len(siCat) - len(np.where(rowsWithGoodSersic)[0])
        self.log.info("Removing %d rows of galaxies with nBulge or nDisk outside of %0.2f <= n <= %0.2f",
                      numRowsNotUsed, minN, maxN)
        siCat = siCat[rowsWithGoodSersic]

        if self.config.doSubSelectSources:
            numRowsNotUsed = len(siCat) - len(siCat["select"])
            self.log.info("Removing %d rows which were not designated as template sources", numRowsNotUsed)
            siCat = siCat[siCat["select"]]

        return siCat
