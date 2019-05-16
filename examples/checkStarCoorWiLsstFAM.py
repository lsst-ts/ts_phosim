import os
import numpy as np

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getModulePath


def precondition(phosimDir):

    # Survey information
    instName = "lsst"
    filterType = FilterType.REF
    ra = 20
    decl = 30
    rotSkyPos = 10
    mjd = 59580.0

    # Declare the SkySim()
    skySim = SkySim()

    # Set the observation information
    skySim.setObservationMetaData(ra, decl, rotSkyPos, mjd)

    # Add the interested stars
    sensorName = "R22_S11"
    starId = [0, 1]
    xInpixelInCam = [3200, 400]
    yInPixelInCam = [3800, 700]
    starMag = [15, 16]
    for ii in range(len(starId)):
        skySim.addStarByChipPos(sensorName, starId[ii], xInpixelInCam[ii],
                                yInPixelInCam[ii], starMag[ii])

    # Output the sky information
    outputFilePath = os.path.join(getModulePath(), "output",
                                  "skyLsstFamInfo.txt")
    skySim.exportSkyToFile(outputFilePath)

    # Set the Telescope facade class
    configFilePath = os.path.join(getModulePath(), "configData",
                                  "telescopeConfig", "GT.inst")
    tele = TeleFacade(configFilePath=configFilePath)
    tele.setSubSysConfigDir(phosimDir=phosimDir)
    tele.setSurveyParam(filterType=filterType, boresight=(ra, decl),
                        rotAngInDeg=rotSkyPos, mjd=mjd)
    tele.setInstName(instName)

    return tele, skySim


def main(phosimDir):

    # Set the output directory
    outputDir = os.path.join(getModulePath(), "output")
    outputImgDir = os.path.join(outputDir, "img")
    cmdSettingFile = os.path.join(getModulePath(), "configData", "cmdFile",
                                  "starDefault.cmd")
    instSettingFile = os.path.join(getModulePath(), "configData", "instFile",
                                   "starSingleExp.inst")

    # Get the objects of TeleFacade and SkySim classes
    tele, skySim = precondition(phosimDir)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile,
                                    cmdFileName="star.cmd")

    # Set the intra- and extra-focal related information
    obsIdList = {"-1": 9005000, "1": 9005001}
    instFileNameList = {"-1": "starExtra.inst", "1": "starIntra.inst"}
    logFileNameList = {"-1": "starExtraPhoSim.log", "1": "starIntraPhoSim.log"}

    outputImgDirIntra = os.path.join(outputImgDir, "intra")
    outputImgDirExtra = os.path.join(outputImgDir, "extra")
    outputImgDirList = {"-1": outputImgDirExtra, "1": outputImgDirIntra}

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
            outputDir, skySim, instSettingFile=instSettingFile,
            instFileName=instFileNameList[str(ii)])

        # Get the argument to run the phosim
        logFilePath = os.path.join(outputImgDirList[str(ii)],
                                   logFileNameList[str(ii)])
        argString = tele.getPhoSimArgs(
            instFilePath, extraCommandFile=cmdFilePath, numPro=1,
            outputDir=outputImgDirList[str(ii)], e2ADC=0,
            logFilePath=logFilePath)
        argStringList.append(argString)

    for ii in range(2):
        tele.runPhoSim(argStringList[ii])


if __name__ == "__main__":

    phosimDir = os.path.join(os.sep, "home", "ttsai", "Document", "bitbucket",
                             "phosim_syseng4")
    main(phosimDir)
