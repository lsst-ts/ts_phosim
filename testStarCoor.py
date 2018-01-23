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
    instSettingFile = "./data/instFile/starDefault.inst"

    folderPath2FocalPlane = os.path.join(phosimDir, "data", "lsst")

    # Set the settings
    obsId = 9006000
    aFilter = "g"
    wavelengthInNm = 500
    mjdTime = 59580.0

    # Generate the camera object
    camera = LsstSimMapper().camera

    # Set the ObservationMetaData
    RA = 0
    Dec = 0
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
    dofInUm = np.zeros(50)

    # Camera dx
    dofInUm[6] = 1000
    tele.setDofInUm(dofInUm)

    # Add the star
    sensorName = "R22_S11"
    starId = [0, 1]
    starMag = [17, 17]
    xInpixelInCam = [1000, 2000]
    yInPixelInCam = [1000, 2036]
    for ii in range(len(starId)):
        skySim.addStarByChipPos(camera, obs, sensorName, starId[ii], xInpixelInCam[ii], yInPixelInCam[ii], 
                                starMag[ii], folderPath2FocalPlane)

    # Write the star physical command file
    cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, cmdFileName="star.cmd")

    # Write the star instance file
    # Use the camRot=0 temporally. Need to update it latter.
    instFilePath = tele.writeStarInstFile(outputDir, skySim, obsId, aFilter, boresight=(RA, Dec), 
                                            camRot=cameraRotation, mjd=mjdTime, sedName="sed_500.txt", 
                                            sciSensorOn=True, instSettingFile=instSettingFile, 
                                            instFileName="star.inst")

    # Write the accumulated DOF file
    tele.writeAccDofFile(outputDir)

    # Get the argument to run the phosim
    logFilePath = os.path.join(outputImgDir, "phosimStar.log")
    argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgDir, 
                                    e2ADC=1, logFilePath=logFilePath)
    
    # Run the phosim
    tele.runPhoSim(argString)
