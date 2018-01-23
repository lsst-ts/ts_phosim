import os
import numpy as np

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
    # There is no camera correction for the comcam at this moment
    M1M3dataDir = "./data/M1M3"
    M2dataDir = "./data/M2"

    outputDir = "./output"
    outputImgDir = "./outputImg"
    cmdSettingFile = "./data/cmdFile/starDefault.cmd"

    # Use the ComCam instance setting file
    instSettingFile = "./data/instFile/starComCamDefault.inst"

    configFilePath = "./data/telescopeConfig/GT.inst"

    # Set the settings
    iSim = 6

    # Assume the observation ID does not change to simplify the demonstration
    obsId = 9006000

    aFilter = "g"
    wavelengthInNm = 500
    zAngleInDeg = 27.0912
    boresight = (0, 0)
    mjdTime = 59552.3

    # Instantiate the subsystems
    M1M3 = M1M3Sim()
    M2 = M2Sim()
    phoSimCommu = PhosimCommu()
    skySim = SkySim()

    # Instantiate the telescope
    # There is no camera correction for the comcam at this moment
    tele = TeleFacade(M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # Set the configuration file path
    tele.setConfigFile(configFilePath)

    # Set the subsystem directory
    tele.setSubSysConfigFile(M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir, phosimDir=phosimDir)

    # Set the comcam camera
    tele.setInstName("comcam10")

    # Add the star
    skySim.addStarByRaDecInDeg(0, 0, 0, 17)
    skySim.addStarByRaDecInDeg(1, 0.0033, 0.0033, 18)

    # Set the degree of freedom for intra- and extra-focal images
    # Double check the direction with Bo
    instFileNameList = {"-1": "starExtra.inst", "1": "starIntra.inst"}
    logFileNameList = {"-1": "starExtraPhoSim.log", "1": "starIntraPhoSim.log"}
    outputImgDirIntra = os.path.join(outputImgDir, "Intra")
    outputImgDirExtra = os.path.join(outputImgDir, "Extra")
    outputImgDirList = {"-1": outputImgDirExtra, "1":outputImgDirIntra}
    argStringList = []
    for ii in (-1, 1):

        # Update the telescope degree of freedom
        dofInUm = np.zeros(50)

        # Camera piston (Change the unit from mm to um)
        dofInUm[5] = ii*tele.defocalDisInMm*1e3

        # Set the degree fo freedom on telescope
        tele.setDofInUm(dofInUm)

        # Write the perturbation command file
        # The perturbation condition is the same for two exposures.
        if (ii == -1):
            pertCmdFilePath = tele.writePertBaseOnConfigFile(outputDir, zAngleInDeg=zAngleInDeg, 
                                                                seedNum=iSim, saveResMapFig=True)

            # Write the star physical command file
            cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, 
                                pertFilePath=pertCmdFilePath, cmdFileName="starComCam.cmd")

        # Write the star instance file
        # Set the camRot=0 for the simplification. But actually, there is no correction for ComCam.
        # However, the coordinate is affected by this.
        instFilePath = tele.writeStarInstFile(outputDir, skySim, obsId, aFilter, boresight=boresight, 
                                                camRot=0, mjd=mjdTime, sedName="sed_500.txt", sciSensorOn=True, 
                                                instSettingFile=instSettingFile, instFileName=instFileNameList[str(ii)])

        # Get the argument to run the phosim
        logFilePath = os.path.join(outputImgDirList[str(ii)], logFileNameList[str(ii)])
        argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, 
                                        outputDir=outputImgDirList[str(ii)], e2ADC=0, logFilePath=logFilePath)
        argStringList.append(argString)

    # Run the phosim
    for ii in range(2):
        tele.runPhoSim(argStringList[ii])
