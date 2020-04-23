#!/usr/bin/env python

import os
import argparse

from lsst.ts.wep.Utility import FilterType
from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath
from lsst.ts.phosim.PlotUtil import plotFwhmOfIters


def main(phosimDir, numPro, iterNum, baseOutputDir, rotCamInDeg=0.0):

    # Survey parameters
    filterType = FilterType.REF

    # Prepare the components
    phosimCmpt = _preparePhosimCmpt(phosimDir, filterType, 0.0, numPro)
    ofcCalc = _prepareOfcCalc(filterType, rotCamInDeg)

    # Set the telescope state to be the same as the OFC
    state0 = ofcCalc.getStateAggregated()
    phosimCmpt.setDofInUm(state0)

    # Do the iteration
    obsId = 9006000
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
        # Rotate OPD in the reversed direction of camera
        phosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
                                        rotOpdInDeg=-rotCamInDeg,
                                        pssnFileName=opdPssnFileName)

        # Get the PSSN from file
        pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
        print("Calculated PSSN is %s." % pssn)

        # Get the GQ effective FWHM from file
        gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
        print("GQ effective FWHM is %.4f." % gqEffFwhm)

        # Set the FWHM data
        refSensorNameList = _getComCamSensorNameList()
        listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
            opdPssnFileName, refSensorNameList)
        ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)

        # Simulate to get the wavefront sensor data from WEP and calculate
        # the DOF
        listOfWfErr = phosimCmpt.mapOpdDataToListOfWfErr(
            opdZkFileName, refSensorNameList)
        ofcCalc.calculateCorrections(listOfWfErr)

        # Set the new aggregated DOF to phosimCmpt
        dofInUm = ofcCalc.getStateAggregated()
        phosimCmpt.setDofInUm(dofInUm)

        # Save the DOF file
        phosimCmpt.saveDofInUmFileForNextIter(
            dofInUm, dofInUmFileName=dofInUmFileName)

        # Add the observation ID by 1
        obsId += 1

    # Summarize the FWHM
    pssnFiles = [os.path.join(baseOutputDir, "%s%d" % (iterDefaultDirName, num),
                 outputImgDirName, opdPssnFileName) for num in range(iterNum)]
    saveToFilePath = os.path.join(baseOutputDir, "fwhmIters.png")
    plotFwhmOfIters(pssnFiles, saveToFilePath=saveToFilePath)


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

    # Update the number of processor if necessary
    if (numPro > 1):
        settingFile = phosimCmpt.getSettingFile()
        settingFile.updateSetting("numPro", numPro)

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

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop simulation in OPD level.")
    parser.add_argument("--numOfProc", type=int, default=1,
                        help="number of processor to run PhoSim (default: 1)")
    parser.add_argument("--iterNum", type=int, default=5,
                        help="number of closed-loop iteration (default: 5)")
    parser.add_argument("--output", type=str, default="",
                        help="output directory")
    parser.add_argument("--rotCam", type=float, default=0.0,
                        help="Rotate camera (degree) in counter-clockwise direction (default: 0.0)")
    args = parser.parse_args()

    # Run the simulation
    phosimDir = getPhoSimPath()

    if (args.output == ""):
        outputDir = getAoclcOutputPath()
    else:
        outputDir = args.output
    os.makedirs(outputDir, exist_ok=True)

    main(phosimDir, args.numOfProc, args.iterNum, outputDir,
         rotCamInDeg=args.rotCam)
