import os
import numpy as np

from wepPhoSim.CamSim import CamSim
from wepPhoSim.M2Sim import M2Sim
from wepPhoSim.M1M3Sim import M1M3Sim
from wepPhoSim.TeleFacade import TeleFacade
from wepPhoSim.PhosimCommu import PhosimCommu

from wepPhoSim.SkySim import SkySim

def main():
    pass

if __name__ == "__main__":

    # PhoSim directory
    phosimDir = "/home/ttsai/Document/bitbucket/phosim_syseng2"

    # Subsystem data direction
    camDataDir = "./data/camera"
    M1M3dataDir = "./data/M1M3"
    M2dataDir = "./data/M2"

    outputDir = "./output"
    outputImgDir = "./outputImg"
    cmdSettingFile = "./data/cmdFile/starDefault.cmd"
    instSettingFile = "./data/instFile/starDefault.inst"

    configFilePath = "./data/telescopeConfig/GT.inst"

    # Set the settings
    iSim = 6
    obsId = 9006000
    aFilter = "g"
    wavelengthInNm = 500
    zAngleInDeg = 27.0912
    rotAngInDeg = -1.2323/np.pi*180.0
    boresight = (0, 0)
    mjdTime = 59552.3

    # Instantiate the subsystems
    cam = CamSim()
    M1M3 = M1M3Sim()
    M2 = M2Sim()
    phoSimCommu = PhosimCommu()
    skySim = SkySim()

    # Instantiate the telescope
    tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # Set the configuration file path
    tele.setConfigFile(configFilePath)

    # Set the subsystem directory
    tele.setSubSysConfigFile(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir, 
                             phosimDir=phosimDir)

    # Update the telescope degree of freedom
    # dofInUm = np.zeros(50)

    # Camera dx
    # dofInUm[6] = 1000
    # tele.setDofInUm(dofInUm)

    # Write the perturbation command file
    pertCmdFilePath = tele.writePertBaseOnConfigFile(outputDir, zAngleInDeg=zAngleInDeg, 
                                rotAngInDeg=rotAngInDeg, seedNum=iSim, saveResMapFig=True)

    # Add the star
    skySim.addStarByRaDecInDeg(0, 1.196, 1.176, 17)
    skySim.addStarByRaDecInDeg(1, 1.196, 1.1793, 16)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, 
                                    pertFilePath=pertCmdFilePath, cmdFileName="star.cmd")

    # Write the star instance file
    # Use the camRot=0 temporally. Need to update it latter.
    instFilePath = tele.writeStarInstFile(outputDir, skySim, obsId, aFilter, boresight=boresight, 
                                            camRot=0, mjd=mjdTime, sedName="sed_500.txt", wfSensorOn=True, 
                                            instSettingFile=instSettingFile, instFileName="star.inst")

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Get the argument to run the phosim
    logFilePath = os.path.join(outputImgDir, "phosimStar.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgDir, 
                                    e2ADC=0, logFilePath=logFilePath)
    
    # Run the phosim
    tele.runPhoSim(argString)
