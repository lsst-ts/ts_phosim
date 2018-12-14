import os
import numpy as np

from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.TeleFacade import TeleFacade
from lsst.ts.phosim.Utility import getModulePath, FilterType


def main(phosimDir):

    # Settings
    outputDir = os.path.join(getModulePath(), "output")
    outputImgDir = os.path.join(outputDir, "img")

    cmdSettingFile = os.path.join(getModulePath(), "configData", "cmdFile",
                                  "opdDefault.cmd")
    instSettingFile = os.path.join(getModulePath(), "configData", "instFile",
                                   "opdDefault.inst")

    # Declare the opd metrology and add the interested field points
    metr = OpdMetrology()
    metr.addFieldXYbyDeg(0, 0)
    metr.addFieldXYbyDeg(0.2, 0.3)

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

    # Set the survey parameters
    obsId = 9006000
    filterType = FilterType.G
    zAngleInDeg = 27.0912
    rotAngInDeg = np.rad2deg(-1.2323)
    tele.setSurveyParam(obsId=obsId, filterType=filterType,
                        zAngleInDeg=zAngleInDeg, rotAngInDeg=rotAngInDeg)

    # Generate the perturbation
    iSim = 6
    pertCmdFilePath = tele.writePertBaseOnConfigFile(outputDir, seedNum=iSim,
                                                     saveResMapFig=True)

    # Update the telescope degree of freedom
    dofInUm = np.zeros(50)

    # Camera dx
    dofInUm[6] = 1000
    tele.accDofInUm(dofInUm)

    # Write the physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile,
                                    pertFilePath=pertCmdFilePath,
                                    cmdFileName="opd.cmd")

    # Write the instance file
    instFilePath = tele.writeOpdInstFile(outputDir, metr,
                                         instSettingFile=instSettingFile,
                                         instFileName="opd.inst")

    # Get the argument to run the PhoSim
    logFilePath = os.path.join(outputImgDir, "opdPhoSim.log")
    argString = tele.getPhoSimArgs(instFilePath, extraCommandFile=cmdFilePath,
                                   numPro=2, outputDir=outputImgDir, e2ADC=0,
                                   logFilePath=logFilePath)

    # Run the PhoSim
    tele.runPhoSim(argString)

    # Analyze the OPD fits images
    opdFitsFile = os.path.join(outputImgDir, "opd_9006000_0.fits.gz")
    zk = metr.getZkFromOpd(opdFitsFile=opdFitsFile)[0]
    print(zk)

    wavelengthInUm = tele.WAVELENGTH_IN_NM * 1e-3
    pssn = metr.calcPSSN(wavelengthInUm, opdFitsFile=opdFitsFile)

    print(pssn)


if __name__ == "__main__":

    phosimDir = os.path.join(os.sep, "home", "ttsai", "Document", "bitbucket",
                             "phosim_syseng4")
    main(phosimDir)
