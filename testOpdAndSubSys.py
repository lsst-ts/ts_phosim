import os
import numpy as np

from wepPhoSim.CamSim import CamSim
from wepPhoSim.M2Sim import M2Sim
from wepPhoSim.M1M3Sim import M1M3Sim
from wepPhoSim.TeleFacade import TeleFacade
from wepPhoSim.PhosimCommu import PhosimCommu

from wepPhoSim.OpdMetrology import OpdMetrology

def main():
    pass

if __name__ == "__main__":

    # PhoSim directory
    phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"

    # Subsystem data direction
    camDataDir = "./data/camera"
    M1M3dataDir = "./data/M1M3"
    M2dataDir = "./data/M2"

    outputDir = "./output"
    outputImgDir = "./outputImg"
    cmdSettingFile = "./data/cmdFile/opdDefault.cmd"
    instSettingFile = "./data/instFile/opdDefault.inst"

    configFilePath = "./data/telescopeConfig/GT.inst"

    # Set the settings
    iSim = 6
    obsId = 9006000
    aFilter = "g"
    wavelengthInNm = 500
    zAngleInDeg = 27.0912
    rotAngInDeg = -1.2323/np.pi*180.0

    # Instantiate the subsystems
    cam = CamSim()
    M1M3 = M1M3Sim()
    M2 = M2Sim()
    phoSimCommu = PhosimCommu()
    metr = OpdMetrology()

    # Instantiate the telescope
    tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # Set the configuration file path
    tele.setConfigFile(configFilePath)

    # Set the subsystem directory
    tele.setSubSysConfigFile(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir, 
                             phosimDir=phosimDir)

    # Write the perturbation command file
    pertCmdFilePath = tele.writePertBaseOnConfigFile(outputDir, zAngleInDeg=zAngleInDeg, 
                                rotAngInDeg=rotAngInDeg, seedNum=iSim, saveResMapFig=True)

    # Set the default LSST GQ metrology
    metr.setDefaultLsstGQ()

    # Add the wavefront sensor
    fieldWFSx, fieldWFSy = metr.getDefaultLsstWfsGQ()
    metr.addFieldXYbyDeg(fieldWFSx, fieldWFSy)

    # Write the opd physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, 
                                    pertFilePath=pertCmdFilePath, cmdFileName="opd.cmd")

    # Write the opd instance file
    instFilePath = tele.writeOpdInstFile(outputDir, metr, obsId, aFilter, wavelengthInNm, 
                                    instSettingFile=instSettingFile, instFileName="opd.inst")

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Get the argument to run the phosim
    logFilePath = os.path.join(outputImgDir, "phosimOpd.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgDir, 
                                    e2ADC=0, logFilePath=logFilePath)
    
    # Run the phosim
    tele.runPhoSim(argString)

    # Calculate the PSSN
    opdFitsFile = os.path.join(outputImgDir, "opd_9006000_0.fits.gz")
    pssn = metr.calcPSSN(wavelengthInNm*1e-3, opdFitsFile=opdFitsFile, zen=0, debugLevel=0)
    print(pssn)
    