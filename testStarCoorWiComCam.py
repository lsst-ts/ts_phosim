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

    # Subsystem data direction
    # There is no camera correction for the comcam at this moment

    outputDir = "./output"
    outputImgDir = "./outputImg"
    cmdSettingFile = "./data/cmdFile/starDefault.cmd"

    # Use the ComCam instance setting file
    instSettingFile = "./data/instFile/starSingleExp.inst"

    # Check to use the comcam instead of lsst in the latter time
    folderPath2FocalPlane = os.path.join(phosimDir, "data", "lsst")

    # Sky information
    outputFilePath = "./output/skyComCamInfo.txt"

    # Set the settings

    # Assume the observation ID does not change to simplify the demonstration
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
    obs = ObservationMetaData(pointingRA=RA, pointingDec=Dec, rotSkyPos=cameraRotation, 
                                mjd=mjdTime)

    # Instantiate the subsystems
    phoSimCommu = PhosimCommu()
    skySim = SkySim()

    # Instantiate the telescope
    # There is no camera correction for the comcam at this moment
    tele = TeleFacade(phoSimCommu=phoSimCommu)

    # Set the subsystem directory
    tele.setSubSysConfigFile(phosimDir=phosimDir)

    # Set the comcam camera
    tele.setInstName("comcam10")

    # Add the star
    sensorName = ["R22_S11", "R22_S10"]
    starMag = [15, 15]
    xInpixelInCam = [500, 800]
    yInPixelInCam = [1000, 1300]

    starId = 0
    for sensor in sensorName:
        for ii in range(len(xInpixelInCam)):
            skySim.addStarByChipPos(camera, obs, sensor, starId, xInpixelInCam[ii], yInPixelInCam[ii], 
                                    starMag[ii], folderPath2FocalPlane)
            starId += 1

    # Export sky information
    skySim.exportSkyToFile(outputFilePath)

    # Set the degree of freedom for intra- and extra-focal images
    # Double check the direction with Bo
    obsIdList = {"-1": 9007000, "1": 9007001}
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

            # Write the star physical command file
            cmdFilePath = tele.writeCmdFile(outputDir, cmdSettingFile=cmdSettingFile, 
                                            cmdFileName="starComCam.cmd")

        # Write the star instance file
        # Set the rot=0 for the simplification. But actually, there is no correction for ComCam.
        # However, the coordinate is affected by this.
        instFilePath = tele.writeStarInstFile(outputDir, skySim, obsIdList[str(ii)], aFilter, boresight=(RA, Dec), 
                                                rot=cameraRotation, mjd=mjdTime, sedName="sed_500.txt", sciSensorOn=True, 
                                                instSettingFile=instSettingFile, instFileName=instFileNameList[str(ii)])

        # Get the argument to run the phosim
        logFilePath = os.path.join(outputImgDirList[str(ii)], logFileNameList[str(ii)])
        argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, 
                                        outputDir=outputImgDirList[str(ii)], e2ADC=0, logFilePath=logFilePath)
        argStringList.append(argString)

    # Run the phosim
    for ii in range(2):
        tele.runPhoSim(argStringList[ii])
