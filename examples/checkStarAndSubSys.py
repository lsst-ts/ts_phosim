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

from lsst.ts.wep.Utility import FilterType, CamType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getConfigDir, getPhoSimPath, getAoclcOutputPath


def main(phosimDir):
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
    zAngleInDeg = 27.0912
    rotSkyPos = np.rad2deg(-1.2323)

    # Set the Telescope facade class
    tele = TeleFacade()
    tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
    tele.setPhoSimDir(phosimDir)
    tele.setSurveyParam(
        obsId=obsId,
        filterType=filterType,
        boresight=(ra, decl),
        zAngleInDeg=zAngleInDeg,
        rotAngInDeg=rotSkyPos,
    )
    tele.setInstName(CamType.LsstFamCam)

    # Declare the SkySim()
    skySim = SkySim()

    # Set the observation information
    skySim.setObservationMetaData(ra, decl, rotSkyPos)

    # Add the interested stars
    sensorName = "R22_S11"
    starId = [0, 1]
    xInpixelInCam = [3200, 300]
    yInPixelInCam = [3800, 1000]
    starMag = [15, 12]
    for ii in range(len(starId)):
        skySim.addStarByChipPos(
            sensorName, starId[ii], xInpixelInCam[ii], yInPixelInCam[ii], starMag[ii]
        )

    # Generate the perturbation
    iSim = 6
    pertCmdFilePath = tele.writePertBaseOnConfigFile(
        outputDir, seedNum=iSim, saveResMapFig=True
    )

    # Update the telescope degree of freedom with camera piston in um
    dofInUm = np.zeros(50)
    dofInUm[5] = 1000
    tele.accDofInUm(dofInUm)

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(
        outputDir,
        cmdSettingFile=cmdSettingFile,
        pertFilePath=pertCmdFilePath,
        cmdFileName="star.cmd",
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
        numPro=1,
        outputDir=outputImgDir,
        e2ADC=0,
        logFilePath=logFilePath,
    )

    # Run the PhoSim
    tele.runPhoSim(argString)


if __name__ == "__main__":
    phosimDir = getPhoSimPath()
    main(phosimDir)
