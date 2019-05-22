import os
import numpy as np

from lsst.ts.wep.Utility import FilterType
from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt
from lsst.ts.phosim.Utility import getConfigDir, getPhoSimPath, \
    getAoclcOutputPath


def main(phosimDir, numPro, iterNum):

    # Survey parameters
    filterType = FilterType.REF
    rotAngInDeg = 0.0

    # Prepare the components
    phosimCmpt = _preparePhosimCmpt(phosimDir, filterType, rotAngInDeg, numPro)
    ofcCalc = _prepareOfcCalc(filterType, rotAngInDeg)

    # Set the telescope state to be the same as the OFC
    state0 = ofcCalc.getStateAggregated()
    phosimCmpt.setDofInUm(state0)

    # Get the sensor name of ComCam
    sensorNameList = _getComCamSensorNameList()

    # Do the iteration
    obsId = 9006000
    baseOutputDir = getAoclcOutputPath()
    opdZkFileName = "opd.zer"
    opdPssnFileName = "PSSN.txt"
    outputDirName = "pert"
    outputImgDirName = "img"
    iterDefaultDirName = "iter"
    dofInUmFileName = "dofPertInNextIter.mat"
    for iterCount in range(iterNum):

        # Set the observation Id
        phosimCmpt.setSurveyParam(obsId=obsId)

        # The iteration directory
        iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

        # Set the output directory
        outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
        phosimCmpt.setOutputDir(outputDir)

        # Set the output image directory
        outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                    outputImgDirName)
        phosimCmpt.setOutputImgDir(outputImgDir)

        # Generate the OPD image
        argString = phosimCmpt.getComCamOpdArgsAndFilesForPhoSim()
        phosimCmpt.runPhoSim(argString)

        # Analyze the OPD data
        phosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
                                        pssnFileName=opdPssnFileName)

        # Get the PSSN from file
        pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
        print("Calculated PSSN is %s." % pssn)

        ######################################
        # Update the following part
        ######################################

        # Set the gain value in ofcCalc by pssn
        ofcCalc.setGainByPSSN(pssn, sensorNameList)

        # Get the GQ effective FWHM from file
        gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
        print(gqEffFwhm)

        # Get the OPD zk from file
        opdZkData = phosimCmpt.getZkFromFile(opdZkFileName)

        # Calculate the new DOF by OFC component
        dofInUm = ofcCalc.calcAggDofForPhoSim(opdZkData, sensorNameList)

        # Set the new DOF to phosimCmpt
        phosimCmpt.setDofInUm(dofInUm)

        # Save the DOF file
        phosimCmpt.saveDofInUmFileForNextIter(
                                    dofInUm, dofInUmFileName=dofInUmFileName)

        # Add the observation ID by 1
        obsId += 1


def _preparePhosimCmpt(phosimDir, filterType, rotAngInDeg, numPro):

    # Set the Telescope facade class
    tele = TeleFacade()
    tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
    tele.setPhoSimDir(phosimDir)

    # Prepare the phosim component
    phosimCmpt = PhosimCmpt(tele)

    # Set the telescope survey parameters
    boresight = (0, 0)
    zAngleInDeg = 27.0912
    phosimCmpt.setSurveyParam(filterType=filterType, boresight=boresight,
                              zAngleInDeg=zAngleInDeg, rotAngInDeg=rotAngInDeg)

    # Set the PhoSim parameters
    phosimCmpt.setPhosimParam(numPro, 1)

    # Set the seed number for M1M3 surface
    seedNum = 6
    phosimCmpt.setSeedNum(seedNum)

    return phosimCmpt


def _prepareOfcCalc(filterType, rotAngInDeg):

    ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)

    ofcCalc.setFilter(filterType)
    ofcCalc.setRotAng(rotAngInDeg)
    ofcCalc.setGainByPSSN()

    return ofcCalc


def _getComCamSensorNameList():

    sensorNameList = ["R22_S00", "R22_S01", "R22_S02", "R22_S10", "R22_S11",
                      "R22_S12", "R22_S20", "R22_S21", "R22_S22"]
    return sensorNameList


if __name__ == "__main__":

    # PhoSim directory
    phosimDir = getPhoSimPath()

    # Number of processor
    numPro = 1

    # Iteration number
    iterNum = 5

    main(phosimDir, numPro, iterNum)
