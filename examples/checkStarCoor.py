import os
import numpy as np

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.TeleFacade import TeleFacade
from lsst.ts.phosim.Utility import getModulePath, FilterType


def main(phosimDir):

    # Settings
    outputDir = os.path.join(getModulePath(), "output")
    outputImgDir = os.path.join(outputDir, "img")

    cmdSettingFile = os.path.join(getModulePath(), "configData", "cmdFile",
                                  "starDefault.cmd")
    instSettingFile = os.path.join(getModulePath(), "configData", "instFile",
                                   "starDefault.inst")

    # Survey information
    obsId = 9006000
    instName = "lsst"
    filterType = FilterType.G
    ra = 20
    decl = 30
    rotSkyPos = 10
    mjd = 59580.0

    # Declare the SkySim()
    skySim = SkySim()

    # Set the focal plane information
    folderPath2FocalPlane = os.path.join(phosimDir, "data", instName)
    skySim.setFolderPath2FocalPlane(folderPath2FocalPlane)

    # Set the observation information
    skySim.setObservationMetaData(ra, decl, rotSkyPos, mjd)

    # Add the interested stars
    sensorName = "R22_S11"
    starId = [0]
    xInpixelInCam = [3200]
    yInPixelInCam = [3800]
    starMag = [15]
    for ii in range(len(starId)):
        skySim.addStarByChipPos(sensorName, starId[ii], xInpixelInCam[ii],
                                yInPixelInCam[ii], starMag[ii])

    # Set the Telescope facade class
    configFilePath = os.path.join(getModulePath(), "configData",
                                  "telescopeConfig", "GT.inst")
    tele = TeleFacade(configFilePath=configFilePath)
    tele.setSubSysConfigDir(phosimDir=phosimDir)
    tele.setSurveyParam(obsId=obsId, filterType=filterType,
                        boresight=(ra, decl), rotAngInDeg=rotSkyPos, mjd=mjd)
    tele.setInstName(instName)

    # Update the telescope degree of freedom
    dofInUm = np.zeros(50)

    # Camera piston in um
    dofInUm[5] = 1000

    tele.accDofInUm(dofInUm)

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile,
                                    cmdFileName="star.cmd")

    # Write the instance file
    instFilePath = tele.writeStarInstFile(outputDir, skySim,
                                          instSettingFile=instSettingFile,
                                          instFileName="star.inst")

    # Get the argument to run the PhoSim
    logFilePath = os.path.join(outputImgDir, "phosimStar.log")
    argString = tele.getPhoSimArgs(instFilePath, extraCommandFile=cmdFilePath,
                                   numPro=1, outputDir=outputImgDir, e2ADC=0,
                                   logFilePath=logFilePath)

    # Run the PhoSim
    tele.runPhoSim(argString)


if __name__ == "__main__":

    phosimDir = os.path.join(os.sep, "home", "ttsai", "Document", "bitbucket",
                             "phosim_syseng4")
    main(phosimDir)
