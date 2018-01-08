import os
import numpy as np

from wepPhoSim.OpdMetrology import OpdMetrology
from wepPhoSim.PhosimCommu import PhosimCommu
from wepPhoSim.TeleFacade import TeleFacade

def main():
    pass

if __name__ == "__main__":

    # Settings 
    phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"
    outputDir = "./output"
    outputImgDir = "./outputImg"
    cmdSettingFile = "./data/cmdFile/opdDefault.cmd"
    instSettingFile = "./data/instFile/opdDefault.inst"

    obsId = 9006050
    aFilter = "g"
    opdWaveInNm = 500.0

    # Declare the opd metrology and add the interested field points
    metr = OpdMetrology()
    metr.addFieldXYbyDeg(0, 0)
    metr.addFieldXYbyDeg(0.2, 0.3)

    # Set the Telescope facade class
    phoSimCommu = PhosimCommu()
    tele = TeleFacade(phoSimCommu=phoSimCommu)
    tele.setSubSysConfigFile(phosimDir=phosimDir)

    # Update the telescope degree of freedom
    dofInUm = np.zeros(50)

    # Camera dx
    dofInUm[6] = 1000
    tele.setDofInUm(dofInUm)

    # Write the physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, 
                                    cmdFileName="opd.cmd")

    # Write the instance file
    instFilePath = tele.writeDefaultOpdInstFile(outputDir, metr, obsId, aFilter, opdWaveInNm, 
                                        instSettingFile=instSettingFile, instFileName="opd.inst")

    # Get the argument to run the PhoSim
    logFilePath = os.path.join(outputImgDir, "opdPhoSim.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, 
                            outputDir=outputImgDir, e2ADC=0, logFilePath=logFilePath)

    # Run the PhoSim
    tele.runPhoSim(argString)

    # Analyze the OPD fits images
    opdFitsFile = os.path.join(outputImgDir, "opd_9006050_0.fits.gz")
    zk = metr.getZkFromOpd(opdFitsFile)[0]
    pssn = metr.calcPSSN(0.5, opdFitsFile)

    print(pssn)




