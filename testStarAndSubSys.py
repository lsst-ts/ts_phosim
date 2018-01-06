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

    # Instantiate the subsystems
    cam = CamSim()
    M1M3 = M1M3Sim()
    M2 = M2Sim()
    phoSimCommu = PhosimCommu()
    skySim = SkySim()

    # Instantiate the telescope
    tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # Subsystem data direction
    camDataDir = "./data/camera"
    M1M3dataDir = "./data/M1M3"
    M2dataDir = "./data/M2"
    phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"

    # Set the configuration file path
    configFilePath = "./data/telescopeConfig/GT.inst"
    tele.setConfigFile(configFilePath)

    # Set the subsystem directory
    tele.setSubSysConfigFile(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir, 
                             phosimDir=phosimDir)

    # Set the path of perturbation file
    iSim = 6
    outputFileDir = "./output"

    # Write the perturbation command file
    pertCmdFilePath = tele.writePertBaseOnConfigFile(outputFileDir, seedNum=iSim, saveResMapFig=True)

    # Add the star
    skySim.addStarByRaDecInDeg(0, 1.196, 1.176, 17)

    # Write the default star command file
    cmdFilePath = tele.writeDefaultStarCmdFile(outputFileDir, pertFilePath=pertCmdFilePath)


    # Write the default star instance file
    obsId = 9006000
    aFilter = "g"
    instFilePath = tele.writeDefaultStarInstFile(outputFileDir, skySim, obsId, aFilter, wfSensorOn=True, 
                                                    instFileName="star.inst")

    # Get the argument to run the phosim
    outputImgFileDir = "./outputImg"
    logFilePath = os.path.join(outputImgFileDir, "phosimStar.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgFileDir, 
                                e2ADC=0, logFilePath=logFilePath)

    # Run the phosim
    tele.runPhoSim(argString)

