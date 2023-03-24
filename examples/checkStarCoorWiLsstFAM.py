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
import numpy as np

from lsst.ts.wep.utility import FilterType, CamType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getConfigDir, getPhoSimPath, getAoclcOutputPath


def main(phosimDir):
    # Set the output directory
    outputDir = getAoclcOutputPath()
    outputImgDir = os.path.join(outputDir, "img")
    os.makedirs(outputImgDir, exist_ok=True)

    configDir = getConfigDir()
    cmdSettingFile = os.path.join(configDir, "cmdFile", "starDefault.cmd")
    instSettingFile = os.path.join(configDir, "instFile", "starSingleExp.inst")

    # Get the objects of TeleFacade and SkySim classes
    tele, skySim = _prepareTeleAndSky(phosimDir)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(
        outputDir, cmdSettingFile=cmdSettingFile, cmdFileName="star.cmd"
    )

    # Set the intra- and extra-focal related information
    obsIdList = {"-1": 9005000, "1": 9005001}
    instFileNameList = {"-1": "starExtra.inst", "1": "starIntra.inst"}
    logFileNameList = {"-1": "starExtraPhoSim.log", "1": "starIntraPhoSim.log"}

    outputImgDirIntra = os.path.join(outputImgDir, "intra")
    outputImgDirExtra = os.path.join(outputImgDir, "extra")
    outputImgDirList = {"-1": outputImgDirExtra, "1": outputImgDirIntra}

    os.makedirs(outputImgDirIntra, exist_ok=True)
    os.makedirs(outputImgDirExtra, exist_ok=True)

    argStringList = []
    for ii in (-1, 1):
        tele.setSurveyParam(obsId=obsIdList[str(ii)])

        # Update the telescope degree of freedom
        dofInUm = np.zeros(50)

        # Camera piston (Change the unit from mm to um)
        dofInUm[5] = ii * tele.getDefocalDistInMm() * 1e3

        # Set the degree fo freedom on telescope
        tele.setDofInUm(dofInUm)

        # Write the star instance file
        instFilePath = tele.writeStarInstFile(
            outputDir,
            skySim,
            instSettingFile=instSettingFile,
            instFileName=instFileNameList[str(ii)],
        )

        # Get the argument to run the phosim
        logFilePath = os.path.join(outputImgDirList[str(ii)], logFileNameList[str(ii)])
        argString = tele.getPhoSimArgs(
            instFilePath,
            extraCommandFile=cmdFilePath,
            numPro=1,
            outputDir=outputImgDirList[str(ii)],
            e2ADC=0,
            logFilePath=logFilePath,
        )
        argStringList.append(argString)

    for ii in range(2):
        tele.runPhoSim(argStringList[ii])


def _prepareTeleAndSky(phosimDir):
    # Survey information
    filterType = FilterType.REF
    ra = 20
    decl = 30
    rotSkyPos = 10

    # Set the Telescope facade class
    tele = TeleFacade()
    tele.setPhoSimDir(phosimDir)
    tele.setSurveyParam(
        filterType=filterType, boresight=(ra, decl), rotAngInDeg=rotSkyPos
    )
    tele.setInstName(CamType.LsstFamCam)

    # Declare the SkySim()
    skySim = SkySim()

    # Set the observation information
    skySim.setObservationMetaData(ra, decl, rotSkyPos)

    # Add the interested stars
    sensorName = "R22_S11"
    starId = [0, 1]
    xInpixelInCam = [3200, 400]
    yInPixelInCam = [3800, 700]
    starMag = [15, 16]
    for ii in range(len(starId)):
        skySim.addStarByChipPos(
            sensorName, starId[ii], xInpixelInCam[ii], yInPixelInCam[ii], starMag[ii]
        )

    # Output the sky information
    outputFilePath = os.path.join(getAoclcOutputPath(), "skyLsstFamInfo.txt")
    skySim.exportSkyToFile(outputFilePath)

    return tele, skySim


if __name__ == "__main__":
    phosimDir = getPhoSimPath()
    main(phosimDir)
