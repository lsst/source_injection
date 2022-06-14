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
Inject sources into coadd-level imaging.
"""

import lsst.pipe.base.connectionTypes as cT
from . import inject_base

__all__ = ["CoaddInjectConfig", "CoaddInjectTask"]


class CoaddInjectConnections(inject_base.BaseInjectConnections,
                             dimensions=("skymap", "tract", "patch", "band"),
                             defaultTemplates={"coaddName": "deep", },
                             ):

    inputExposure = cT.Input(
        doc="Exposure to inject synthetic sources into.",
        name="{coaddName}Coadd",
        storageClass="ExposureF",
        dimensions=("skymap", "tract", "patch", "band"),
    )

    outputExposure = cT.Output(
        doc="Injected Exposure.",
        name="{siPrefix}{coaddName}Coadd",
        storageClass="ExposureF",
        dimensions=("skymap", "tract", "patch", "band"),
    )


class CoaddInjectConfig(inject_base.BaseInjectConfig,
                        pipelineConnections=CoaddInjectConnections,
                        ):

    pass


class CoaddInjectTask(inject_base.BaseInjectTask):
    """Class to inject into coadd-level images.
    """

    _DefaultName = "coaddInjectTask"
    ConfigClass = CoaddInjectConfig
