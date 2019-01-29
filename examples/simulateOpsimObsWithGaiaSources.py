import os

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.TeleFacade import TeleFacade
from lsst.ts.phosim.Utility import getModulePath, getOpsimObservation


def main(phosimDir):

    # Settings
    outputDir = os.path.join(getModulePath(), "output")
    outputImgDir = os.path.join(outputDir, "img")

    cmdSettingFile = os.path.join(getModulePath(), "configData", "cmdFile",
                                  "starDefault.cmd")

    # Get the Opsim observation
    tract = 0
    target = 0
    obs = getOpsimObservation(tract=tract, target=target)

    # Instantiate SkySim populated with Gaia sources
    skySim = SkySim.fromObservationWithGaiaSources(obs)

    # Set the focal plane information
    instName = "lsst"
    folderPath2FocalPlane = os.path.join(phosimDir, "data", instName)
    skySim.setFolderPath2FocalPlane(folderPath2FocalPlane)

    # Export sky information
    outputSkyFilePath = os.path.join(getModulePath(), "output",
                                     "skyWfsInfo.txt")
    skySim.exportSkyToFile(outputSkyFilePath)

    # Set the Telescope facade class
    configFilePath = os.path.join(getModulePath(), "configData",
                                  "telescopeConfig", "GT.inst")
    tele = TeleFacade(configFilePath=configFilePath)
    tele.setSubSysConfigDir(phosimDir=phosimDir)
    tele.setObservation(obs)
    tele.setInstName(instName)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile,
                                    cmdFileName="star.cmd")

    # Write the instance file
    instFilePath = tele.writeStarInstFile(outputDir, skySim, 
        instFileName="star.inst")

    # Get the argument to run the PhoSim
    logFilePath = os.path.join(outputImgDir, "phosimStar.log")
    argString = tele.getPhoSimArgs(instFilePath, extraCommandFile=cmdFilePath,
                                   numPro=8, outputDir=outputImgDir, e2ADC=1,
                                   logFilePath=logFilePath)

    # Run the PhoSim
    tele.runPhoSim(argString)


if __name__ == "__main__":

    phosimDir = "/project/activeoptics/firstdonuts/phosim_syseng4"
    main(phosimDir)
