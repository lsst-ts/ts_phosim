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

from lsst.ts.wep.utils import FilterType, CamType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getConfigDir, getPhoSimPath, getAoclcOutputPath


def main(phosimDir, numPro):
    # Settings
    outputDir = getAoclcOutputPath()
    outputImgDir = os.path.join(outputDir, "img")
    os.makedirs(outputImgDir, exist_ok=True)

    configDir = getConfigDir()
    cmdSettingFile = os.path.join(configDir, "cmdFile", "starDefault.cmd")
    instSettingFile = os.path.join(configDir, "instFile", "starSingleExp.inst")

    # Survey information
    obsId = 9006000
    filterType = FilterType.REF
    ra = 20
    decl = 30
    rotSkyPos = 10

    # Set the Telescope facade class
    tele = TeleFacade()
    tele.setPhoSimDir(phosimDir)
    tele.setSurveyParam(
        obsId=obsId, filterType=filterType, boresight=(ra, decl), rotAngInDeg=rotSkyPos
    )
    tele.setInstName(CamType.LsstCam)

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Declare the SkySim()
    skySim = SkySim()

    # Set the observation information
    skySim.setObservationMetaData(ra, decl, rotSkyPos)

    # Add the interested stars
    sensorName = [
        "R44_SW1",
        "R00_SW0",
        "R44_SW0",
        "R00_SW1",
        "R04_SW0",
        "R40_SW1",
        "R04_SW1",
        "R40_SW0",
    ]
    xInpixelInCam = [500, 800]
    yInPixelInCam = [1000, 1300]
    starMag = [15, 15]
    starId = 0
    for sensor in sensorName:
        for ii in range(len(starMag)):
            skySim.addStarByChipPos(
                sensor, starId, xInpixelInCam[ii], yInPixelInCam[ii], starMag[ii]
            )
            starId += 1

    # Export sky information
    outputSkyFilePath = os.path.join(outputDir, "skyWfsInfo.txt")
    skySim.exportSkyToFile(outputSkyFilePath)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(
        outputDir, cmdSettingFile=cmdSettingFile, cmdFileName="star.cmd"
    )

    # Write the instance file
    instFilePath = tele.writeStarInstFile(
        outputDir, skySim, instSettingFile=instSettingFile, instFileName="star.inst"
    )

    # Get the argument to run the PhoSim
    logFilePath = os.path.join(outputImgDir, "phosimStar.log")
    argString = tele.getPhoSimArgs(
        instFilePath,
        extraCommandFile=cmdFilePath,
        numPro=numPro,
        outputDir=outputImgDir,
        e2ADC=0,
        logFilePath=logFilePath,
    )

    # Run the PhoSim
    tele.runPhoSim(argString)


if __name__ == "__main__":
    phosimDir = getPhoSimPath()
    numPro = 1
    main(phosimDir, numPro)
