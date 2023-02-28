#!/usr/bin/env python

# This file is part of ts_phosim.
#
# Developed for the LSST Telescope and Site Systems.
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

import os
import argparse
import tempfile

from lsst.ts.phosim.utils.CreatePhosimDonutTemplates import CreatePhosimDonutTemplates

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate donut templates for AOS using Phosim."
    )
    parser.add_argument(
        "--numOfProc",
        type=int,
        default=1,
        help="Number of processor to run PhoSim. (default: 1)",
    )
    parser.add_argument(
        "--detectorList",
        type=str,
        default="",
        help="""
        Specify detectors.
        By default will generate templates for all detectors.
        (Example Input: "R22_S00 R22_S01 R22_S02")
        """,
    )
    parser.add_argument(
        "--templateSize",
        type=int,
        default=240,
        help="Size of each side of the template in pixels. (default: 240)",
    )
    parser.add_argument(
        "--intraVisitId",
        type=int,
        default=9006002,
        help="Visit ID for phosim intrafocal images. (default: 9006002)",
    )
    parser.add_argument(
        "--extraVisitId",
        type=int,
        default=9006001,
        help="Visit ID for phosim extrafocal images. (default: 9006001)",
    )
    parser.add_argument(
        "--isrConfigFile",
        type=str,
        default="",
        help="""Path to ISR pipeline configuration file.
        By default will use createPhosimDonutTemplateConfig.yaml
        in ts_phosim/policy/donutTemplateData.
        """,
    )
    parser.add_argument(
        "--tempWorkDir",
        type=str,
        default="",
        help="""Location to put temporary work files and directories.
        By default will use python tempfile to create one in the current
        working directory.""",
    )
    args = parser.parse_args()

    if args.tempWorkDir == "":
        tempTestDirectory = tempfile.TemporaryDirectory(dir=os.environ["PWD"])
    else:
        tempTestDirectory = tempfile.TemporaryDirectory(dir=args.tempWorkDir)
    tempWorkDir = tempTestDirectory.name

    # Run tasks
    phosimDonuts = CreatePhosimDonutTemplates(tempWorkDir)
    phosimDonuts.createWorkDirectories()
    detListPhosim = phosimDonuts.createDetectorLists(detectorStr=args.detectorList)
    phosimDonuts.generateDefocalImages(detListPhosim, args.numOfProc)
    phosimDonuts.repackagePhosimImages()
    phosimDonuts.ingestImages()
    phosimDonuts.runISR(args.isrConfigFile)
    phosimDonuts.cutOutIntraExtraTemplates(
        args.templateSize, args.intraVisitId, args.extraVisitId
    )
    phosimDonuts.cleanUpWorkDirs()
