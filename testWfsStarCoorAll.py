import os
import numpy as np

from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.utils import ObservationMetaData

from wepPhoSim.TeleFacade import TeleFacade
from wepPhoSim.PhosimCommu import PhosimCommu
from wepPhoSim.SkySim import SkySim

def main():
    pass

if __name__ == "__main__":

    # PhoSim directory
    phosimDir = "/home/ttsai/Document/bitbucket/phosim_syseng2"

    outputDir = "./output"
    outputImgDir = "./outputImg"
    cmdSettingFile = "./data/cmdFile/starDefault.cmd"
    instSettingFile = "./data/instFile/starSingleExp.inst"

    folderPath2FocalPlane = os.path.join(phosimDir, "data", "lsst")

    # Sky information
    outputFilePath = "./output/skyWfsAllInfo.txt"

    # Set the settings
    obsId = 9001000
    aFilter = "g"
    wavelengthInNm = 500
    mjdTime = 59580.0

    # Generate the camera object
    camera = LsstSimMapper().camera

    # Set the ObservationMetaData
    RA = 0
    Dec = 0

    # The unit of camera rotation angle is in degree
    cameraRotation = 0
    cameraMJD = 59580.0

    obs = ObservationMetaData(pointingRA=RA, pointingDec=Dec, rotSkyPos=cameraRotation, 
                                mjd=mjdTime)

    # Instantiate the subsystems
    phoSimCommu = PhosimCommu()
    skySim = SkySim()

    # Instantiate the telescope
    tele = TeleFacade(phoSimCommu=phoSimCommu)

    # Set the subsystem directory
    tele.setSubSysConfigFile(phosimDir=phosimDir)

    # Update the telescope degree of freedom
    # dofInUm = np.zeros(50)

    # Camera piston in um
    # dofInUm[5] = 1000

    # Camera dx in um
    # dofInUm[6] = 500

    # Set the telescope DOF
    # tele.setDofInUm(dofInUm)

    # Add the star on WFS
    sensorName = ["R44_S00_C0", "R00_S22_C1", "R44_S00_C1", "R00_S22_C0", "R04_S20_C1", 
                    "R40_S02_C0", "R04_S20_C0", "R40_S02_C1"]
    starId = 0
    numOfStar = 20
    starMag = 15*np.ones(numOfStar)
    for sensor in sensorName:

        xInPixelInCam = np.random.rand(numOfStar)*2000
        yInPixelInCam = np.random.rand(numOfStar)*4072

        xInPixelInCam = xInPixelInCam.astype(int)
        yInPixelInCam = yInPixelInCam.astype(int)

        for ii in range(numOfStar):
            skySim.addStarByChipPos(camera, obs, sensor, starId, xInPixelInCam[ii], yInPixelInCam[ii], 
                                    starMag[ii], folderPath2FocalPlane)
            starId += 1

    # Export sky information
    skySim.exportSkyToFile(outputFilePath)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, cmdFileName="starWfsAll.cmd")

    # Write the star instance file
    # Use the rot=0 temporally. Need to update it latter.
    instFilePath = tele.writeStarInstFile(outputDir, skySim, obsId, aFilter, boresight=(RA, Dec), 
                                            rot=cameraRotation, mjd=mjdTime, sedName="sed_500.txt", 
                                            wfSensorOn=True, instSettingFile=instSettingFile, 
                                            instFileName="starWfsAll.inst")

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Get the argument to run the phosim
    logFilePath = os.path.join(outputImgDir, "phosimStarWfsAll.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=8, outputDir=outputImgDir, 
                                    e2ADC=1, logFilePath=logFilePath)
    
    # Run the phosim
    tele.runPhoSim(argString)
