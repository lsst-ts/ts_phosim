import os
import numpy as np

from wepPhoSim.CamSim import CamSim
from wepPhoSim.M2Sim import M2Sim
from wepPhoSim.M1M3Sim import M1M3Sim
from wepPhoSim.TeleFacade import TeleFacade

from wepPhoSim.OpdMetrology import OpdMetrology

def main():
    pass

if __name__ == "__main__":

    cam = CamSim()
    M1M3 = M1M3Sim()
    M2 = M2Sim()
    phoSimCommu = PhosimCommu()
    metr = OpdMetrology()

    tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # Set the telescope degree of freedom

    # Subsystem data direction
    camDataDir = "../data/camera"
    M1M3dataDir = "../data/M1M3"
    M2dataDir = "../data/M2"
    phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"

    # Set the configuration file path
    configFilePath = "../data/telescopeConfig/GT.inst"
    tele.setConfigFile(configFilePath)

    # Get the value
    varName = "zenithAngle"
    value = tele.getConfigValue(varName, index=1)

    # Set the subsystem directory
    tele.setSubSysConfigFile(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir, 
                             phosimDir=phosimDir)

    # Set the path of perturbation file
    iSim = 6
    outputFileDir = "../output"

    # Write the perturbation command file
    # content = tele.writePertBaseOnConfigFile(outputFileDir, seedNum=iSim, saveResMapFig=True)

    # Set the default LSST GQ metrology
    metr.setDefaultLsstGQ()

    # Add the wavefront sensor
    fieldWFSx, fieldWFSy = metr.getDefaultLsstWfsGQ()
    metr.addFieldXYbyDeg(fieldWFSx, fieldWFSy)

    # Write the opd physical command file
    pertFilePath = os.path.join(outputFileDir, "pert.cmd")
    cmdFilePath = tele.writeDefaultOpdCmdFile(outputFileDir, pertFilePath=pertFilePath)

    # Write the opd instance file
    obsId = 9006000
    aFilter = "g"
    wavelengthInNm = 500
    instFilePath = tele.writeDefaultOpdInstFile(outputFileDir, metr, obsId, aFilter, wavelengthInNm)

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputFileDir)

    # Get the argument to run the phosim
    outputImgFileDir = "../outputImg"
    logFilePath = os.path.join(outputImgFileDir, "phosimOpd.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgFileDir, 
                                    e2ADC=0, logFilePath=logFilePath)
    
    # Run the phosim
    # tele.runPhoSim(argString)

    # Calculate the PSSN
    zen = tele.getConfigValue("zenithAngle", index=1)
    opdFitsFile = os.path.join(outputImgFileDir, "opd_9006000_0.fits.gz")
    wavelengthInUm = 0.5
    pssn = metr.calcPSSN(opdFitsFile, wavelengthInUm, zen=zen, debugLevel=0)
    