import os

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.TeleFacade import TeleFacade
from lsst.ts.phosim.Utility import getModulePath


def main(phosimDir):

    # Settings
    outputDir = os.path.join(getModulePath(), "output")
    outputImgDir = os.path.join(outputDir, "img")

    cmdSettingFile = os.path.join(getModulePath(), "configData", "cmdFile",
                                  "starDefault.cmd")
    instSettingFile = os.path.join(getModulePath(), "configData", "instFile",
                                   "starSingleExp.inst")

    # Survey information
    obsId = 9006000
    instName = "lsst"
    filterType = FilterType.REF
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
    sensorName = ["R44_S00_C0", "R00_S22_C1", "R44_S00_C1", "R00_S22_C0",
                  "R04_S20_C1", "R40_S02_C0", "R04_S20_C0", "R40_S02_C1"]
    xInpixelInCam = [500, 800]
    yInPixelInCam = [1000, 1300]
    starMag = [15, 15]
    starId = 0
    for sensor in sensorName:
        for ii in range(len(starMag)):
            skySim.addStarByChipPos(sensor, starId, xInpixelInCam[ii],
                                    yInPixelInCam[ii], starMag[ii])
            starId += 1

    # Export sky information
    outputSkyFilePath = os.path.join(getModulePath(), "output",
                                     "skyWfsInfo.txt")
    skySim.exportSkyToFile(outputSkyFilePath)

    # Set the Telescope facade class
    configFilePath = os.path.join(getModulePath(), "configData",
                                  "telescopeConfig", "GT.inst")
    tele = TeleFacade(configFilePath=configFilePath)
    tele.setSubSysConfigDir(phosimDir=phosimDir)
    tele.setSurveyParam(obsId=obsId, filterType=filterType,
                        boresight=(ra, decl), rotAngInDeg=rotSkyPos, mjd=mjd)
    tele.setInstName(instName)

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
                                   numPro=8, outputDir=outputImgDir, e2ADC=0,
                                   logFilePath=logFilePath)

    # Run the PhoSim
    tele.runPhoSim(argString)


if __name__ == "__main__":

    phosimDir = os.path.join(os.sep, "home", "ttsai", "Document", "bitbucket",
                             "phosim_syseng4")
    main(phosimDir)
