import os

from lsst.ts.wep.Utility import FilterType, CamType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getConfigDir, getPhoSimPath, \
    getAoclcOutputPath


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
    tele.setSurveyParam(obsId=obsId, filterType=filterType,
                        boresight=(ra, decl), rotAngInDeg=rotSkyPos)
    tele.setInstName(CamType.LsstCam)

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Declare the SkySim()
    skySim = SkySim()

    # Set the observation information
    mjd = tele.getCamMjd()
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
    outputSkyFilePath = os.path.join(outputDir, "skyWfsInfo.txt")
    skySim.exportSkyToFile(outputSkyFilePath)

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
                                   numPro=numPro, outputDir=outputImgDir,
                                   e2ADC=0, logFilePath=logFilePath)

    # Run the PhoSim
    tele.runPhoSim(argString)


if __name__ == "__main__":

    phosimDir = getPhoSimPath()
    numPro = 1
    main(phosimDir, numPro)
