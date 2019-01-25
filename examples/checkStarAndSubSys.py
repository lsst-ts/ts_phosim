import os
import numpy as np

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.TeleFacade import TeleFacade
from lsst.ts.phosim.Utility import getModulePath, createObservation


def main(phosimDir):

    # Settings
    outputDir = os.path.join(getModulePath(), "output")
    outputImgDir = os.path.join(outputDir, "img")

    cmdSettingFile = os.path.join(getModulePath(), "configData", "cmdFile",
                                  "starDefault.cmd")
    instOverrideFile = os.path.join(getModulePath(), "configData", "instFile",
                                   "starSingleExp.inst")

    # Survey information
    obsId = 9006000
    instName = "lsst"
    filterType = FilterType.REF
    zAngleInDeg = 27.0912

    ra = 20
    decl = 30
    rotSkyPos = np.rad2deg(-1.2323)
    mjd = 59580.0

    obs = createObservation(obsId=obsId, filterType=filterType,
                        boresight=(ra, decl), zAngleInDeg=zAngleInDeg,
                        rotAngInDeg=rotSkyPos, mjd=mjd)

    # Declare the SkySim()
    skySim = SkySim()

    # Set the focal plane information
    folderPath2FocalPlane = os.path.join(phosimDir, "data", instName)
    skySim.setFolderPath2FocalPlane(folderPath2FocalPlane)

    # Set the observation information
    skySim.setObservationMetaData(obs)

    # Add the interested stars
    sensorName = "R22_S11"
    starId = [0, 1]
    xInpixelInCam = [3200, 300]
    yInPixelInCam = [3800, 1000]
    starMag = [15, 12]
    for ii in range(len(starId)):
        skySim.addStarByChipPos(sensorName, starId[ii], xInpixelInCam[ii],
                                yInPixelInCam[ii], starMag[ii])

    # Set the Telescope facade class
    configFilePath = os.path.join(getModulePath(), "configData",
                                  "telescopeConfig", "GT.inst")
    tele = TeleFacade(configFilePath=configFilePath)

    # Subsystem data direction
    camDataDir = os.path.join(getModulePath(), "configData", "camera")
    M1M3dataDir = os.path.join(getModulePath(), "configData", "M1M3")
    M2dataDir = os.path.join(getModulePath(), "configData", "M2")
    tele.setSubSysConfigDir(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir,
                            M2dataDir=M2dataDir, phosimDir=phosimDir)

    tele.setObservation(obs)
    tele.setInstName(instName)

    # Generate the perturbation
    iSim = 6
    pertCmdFilePath = tele.writePertBaseOnConfigFile(outputDir, seedNum=iSim,
                                                     saveResMapFig=True)

    # Update the telescope degree of freedom
    dofInUm = np.zeros(50)

    # Camera piston in um
    dofInUm[5] = 1000

    tele.accDofInUm(dofInUm)

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile,
                                    pertFilePath=pertCmdFilePath,
                                    cmdFileName="star.cmd")

    # Write the instance file
    instFilePath = tele.writeStarInstFile(outputDir, skySim,
                                          instOverrideFile=instOverrideFile,
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
