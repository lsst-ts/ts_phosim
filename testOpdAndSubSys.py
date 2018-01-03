import os

from wepPhoSim.OpdMetrology import OpdMetrology
from wepPhoSim.PhosimCommu import PhosimCommu

from wepPhoSim.CamSim import CamSim
from wepPhoSim.M2Sim import M2Sim
from wepPhoSim.M1M3Sim import M1M3Sim

def main(phosimDir, outputDir, obsId, aFilter, opdWaveInNm, numProc=1):
    
    # Declare the opd metrology call
    metr = OpdMetrology()

    # Set the lsst camera
    metr.setDefaultLsstGQ()

    # Write the opd instance
    phosimCom = PhosimCommu(phosimDir)
    instFilePath = os.path.join(outputDir, "opd.inst")
    
    aFilterId = phosimCom.getFilterId(aFilter)
    phosimCom.getDefaultOpdInstance(obsId, aFilterId, filePath=instFilePath)

    opdSourceContent = ""
    for ii in range(len(metr.fieldX)):
        opdSourceContent += phosimCom.generateOpd(ii, metr.fieldX[ii], metr.fieldY[ii], opdWaveInNm)
    phosimCom.writeToFile(instFilePath, opdSourceContent, mode="a")

    # Write the opd physical command
    cmdFilePath = os.path.join(outputDir, "opd.cmd")
    phosimCom.getDefaultOpdCmd(filePath=cmdFilePath)

    # Prepare the argument to run the phosim
    logFilePath = os.path.join(outputDir, "opdPhoSim.log")
    argString = phosimCom.getPhoSimArgs(instFilePath, extraCommand=cmdFilePath, numProc=numProc, outputDir=outputDir, 
                                        instrument="lsst", e2ADC=0, logFilePath=logFilePath)

    # Run the PhoSim
    # phosimCom.runPhoSim(argstring=argString)
    print(argString)

if __name__ == "__main__":

    phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"
    outputDir = "./output"
    opdWaveInNm = 500.0

    aFilter = "g"
    obsId = 9006050
    numProc = 2

    zenithAngleInDeg = 27.0912
    zAngleInRad = zenithAngleInDeg/180*np.pi

    # Add camera
    camRotationInRadian = -1.2323
    tempInDegC = 6.5650

    camDataDir = "./data/camera"
    camSim = CamSim()
    camSim.setCamDataDir(self.camDataDir)
    camSim.setRotAngInRad(camRotationInRadian)
    camSim.setBodyTempInDegC(tempInDegC)

    distortionInMn = camSim.getCamDistortionInMm(zAngleInRad, distType)

    main(phosimDir, outputDir, obsId, aFilter, opdWaveInNm, numProc=numProc)