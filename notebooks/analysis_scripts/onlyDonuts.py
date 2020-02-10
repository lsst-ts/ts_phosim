import os
import argparse
import numpy as np
import shutil

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import FilterType, CamType, runProgram
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory
from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, \
                                   getConfigDir
from lsst.ts.phosim.PlotUtil import plotFwhmOfIters

from lsst.ts.wep.bsc.Filter import Filter

from baseComcamLoop import AnalysisPhosimCmpt
from createPhosimCatalog import createPhosimCatalog

class genDonuts():

    def main(self, phosimDir, numPro, baseOutputDir, 
            testName, isEimg=False, genFlats=True,
            surveyFilter=None, starMag=15, 
            useMinDofIdx=False, inputSkyFilePath="", m1m3ForceError=0.05):
            
        # Make the ISR directory
        isrDirName = "input"
        isrDir = os.path.join(baseOutputDir, isrDirName)
        if genFlats is True:
            self._makeDir(isrDir)

        # Survey parameters
        surveySettingFilePath = os.path.join(getConfigDir(),
                                            "surveySettings.yaml")
        surveySettings = ParamReader(filePath=surveySettingFilePath)
        if surveyFilter is None:
            filterType = FilterType.fromString(
                surveySettings.getSetting("filterType"))
        else:
            filterType = FilterType.fromString(surveyFilter)
        raInDeg = surveySettings.getSetting("raInDeg")
        decInDeg = surveySettings.getSetting("decInDeg")
        rotAngInDeg = surveySettings.getSetting("rotAngInDeg")

        # Prepare the components
        phosimCmpt = self._preparePhosimCmpt(phosimDir, filterType, raInDeg, decInDeg,
                                        rotAngInDeg, numPro, isEimg,
                                        m1m3ForceError)

        wepCalc = self._prepareWepCalc(isrDir, filterType, raInDeg, decInDeg,
                                rotAngInDeg, isEimg)

        tele = phosimCmpt.getTele()
        defocalDisInMm = tele.getDefocalDistInMm()
        wepCalc.setDefocalDisInMm(defocalDisInMm)

        ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg)

        # Only use 10 hexapod positions and first 3 bending modes of M1M3 and M2
        if (useMinDofIdx):
            self._useMinDofIdx(ofcCalc)

        # Set the telescope state to be the same as the OFC
        state0 = ofcCalc.getStateAggregated()
        phosimCmpt.setDofInUm(state0)

        # Do the iteration
        filt_num_dict = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}
        obsId = 9006000 + 10*filt_num_dict[surveyFilter]
        opdZkFileName = str("opd.zer" + '.' + testName)
        wfsZkFileName = str("wfs.zer" + '.' + testName)
        opdPssnFileName = "PSSN.txt"
        outputDirName = "pert"
        outputImgDirName = "img"
        iterDefaultDirName = "iter"
        dofInUmFileName = "dofPertInNextIter.mat"
        skyInfoFileName = "skyComCamInfo.txt"

        # Set the observation Id
        phosimCmpt.setSurveyParam(obsId=obsId)

        # The iteration directory
        iterDirName = "%s_%.2f" % (surveyFilter, starMag)
        if os.path.isdir(os.path.join(baseOutputDir, iterDirName)):
            self._eraseFolderContent(os.path.join(baseOutputDir,
                                                  iterDirName))

        # Set the output directory
        outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
        phosimCmpt.setOutputDir(outputDir)

        # Set the output image directory
        outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                    outputImgDirName)
        phosimCmpt.setOutputImgDir(outputImgDir)

        # Prepare the faked sky
        if (inputSkyFilePath == ""):
            # According to the OPD field positions
            metr = phosimCmpt.getOpdMetr()
            skySim = self._prepareSkySim(metr, starMag)
            print("Use the default OPD field positions to be star positions.")
            print("The star magnitude is chosen to be %.2f." % starMag)
        else:
            skySim = self._prepareSkySimBySkyFile(inputSkyFilePath)

        # Output the sky information
        skySim, wepCalc = self._outputSkyInfo(outputDir, skyInfoFileName, skySim, wepCalc)

        # Assign the entra- and intra-focal observation Id
        extraObsId = obsId + 1
        intraObsId = obsId + 2

        # Generate the defocal images
        simSeed = 1000
        argStringList = phosimCmpt.getComCamStarArgsAndFilesForPhoSim(
            extraObsId, intraObsId, skySim, simSeed=simSeed,
            cmdSettingFileName="starDefault.cmd",
            instSettingFileName="starSingleExp.inst")

        for argString in argStringList:
            phosimCmpt.runPhoSim(argString)

        # Repackage the images based on the image type
        if (isEimg):
            phosimCmpt.repackageComCamEimgFromPhoSim()
        else:
            phosimCmpt.repackageComCamAmpImgFromPhoSim()

    def _eraseFolderContent(self, targetDir):

        for theFile in os.listdir(targetDir):
            filePath = os.path.join(targetDir, theFile)
            if os.path.isfile(filePath):
                os.unlink(filePath)
            elif os.path.isdir(filePath):
                shutil.rmtree(filePath)

    def _outputSkyInfo(self, outputDir, skyInfoFileName, skySim, wepCalc):

        outputSkyInfoFilePath = os.path.join(outputDir, skyInfoFileName)
        skySim.exportSkyToFile(outputSkyInfoFilePath)
        wepCalc.setSkyFile(outputSkyInfoFilePath)

        return skySim, wepCalc


    def _preparePhosimCmpt(self, phosimDir, filterType, raInDeg, decInDeg, rotAngInDeg,
                        numPro, isEimg, m1m3ForceError):

        # Set the Telescope facade class
        tele = TeleFacade()
        tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
        tele.setPhoSimDir(phosimDir)

        # Prepare the phosim component
        phosimCmpt = AnalysisPhosimCmpt(tele)

        # Set the telescope survey parameters
        boresight = (raInDeg, decInDeg)
        zAngleInDeg = 27.0912
        phosimCmpt.setSurveyParam(filterType=filterType, boresight=boresight,
                                zAngleInDeg=zAngleInDeg, rotAngInDeg=rotAngInDeg)

        # Update the setting file if needed
        settingFile = phosimCmpt.getSettingFile()
        if (numPro > 1):
            settingFile.updateSetting("numPro", numPro)
        if isEimg:
            settingFile.updateSetting("e2ADC", 0)

        # Set the seed number for M1M3 surface
        seedNum = 6
        phosimCmpt.setSeedNum(seedNum)

        # Set the M1M3 force error
        phosimCmpt.setM1M3ForceError(m1m3ForceError)

        return phosimCmpt


    def _prepareWepCalc(self, isrDirPath, filterType, raInDeg, decInDeg, rotAngInDeg,
                        isEimg):

        wepCalc = WEPCalculationFactory.getCalculator(CamType.ComCam, isrDirPath)
        wepCalc.setFilter(filterType)
        wepCalc.setBoresight(raInDeg, decInDeg)
        wepCalc.setRotAng(rotAngInDeg)

        if (isEimg):
            settingFile = wepCalc.getSettingFile()
            settingFile.updateSetting("imageType", "eimage")

        return wepCalc


    def _prepareOfcCalc(self, filterType, rotAngInDeg):

        ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
        ofcCalc.setFilter(filterType)
        ofcCalc.setRotAng(rotAngInDeg)
        ofcCalc.setGainByPSSN()

        return ofcCalc


    def _prepareSkySim(self, opdMetr, starMag):

        skySim = SkySim()

        starId = 0
        raInDegList, declInDegList = opdMetr.getFieldXY()
        for raInDeg, declInDeg in zip(raInDegList, declInDegList):
            # It is noted that the field position might be < 0. But it is not the
            # same case for ra (0 <= ra <= 360).
            if (raInDeg < 0):
                raInDeg += 360.0
            skySim.addStarByRaDecInDeg(starId, raInDeg, declInDeg, starMag)
            starId += 1

        return skySim


    def _prepareSkySimBySkyFile(self, inputSkyFilePath):

        skySim = SkySim()

        absSkyFilePath = os.path.abspath(inputSkyFilePath)
        skySim.addStarByFile(absSkyFilePath)

        return skySim


    def _useMinDofIdx(self, ofcCalc):

        ztaac = ofcCalc.getZtaac()

        m1m3Bend = np.zeros(20, dtype=int)
        m1m3Bend[0: 3] = 1

        m2Bend = np.zeros(20, dtype=int)
        m2Bend[0: 3] = 1

        ztaac.setZkAndDofInGroups(m1m3Bend=m1m3Bend, m2Bend=m2Bend)

    def _makeDir(self, directory):

        if (not os.path.exists(directory)):
            os.makedirs(directory)

if __name__ == "__main__":

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop simulation (default is amp files).")
    parser.add_argument("--numOfProc", type=int, default=8,
                        help="number of processor to run PhoSim (default: 8)")
    parser.add_argument("--output", type=str, default="",
                        help="output directory")
    parser.add_argument("--testOutput", type=str, default="", help="test output directory")
    parser.add_argument("--testLabel", type=str, default="1",
                        help="filename identifier for test files")
    parser.add_argument('--eimage', default=False, action='store_true',
                        help='Use the eimage files')
    parser.add_argument('--minDof', default=False, action='store_true',
                        help='Use 10 hexapod positions and first 3 bending modes of M1M3 and M2')
    parser.add_argument("--skyFile", type=str, default="",
                        help="Star Id, ra, dec, and magnitude")
    parser.add_argument("--m1m3FErr", type=float, default=0.05,
                        help="Ratio of M1M3 actuator force error between 0 and 1 (default: 0.05)")
    parser.add_argument('--clobber', default=False, action='store_true',
                        help='Delete existing output directory')
    parser.add_argument("--raShift", type=float, default=0.0)
    parser.add_argument("--decShift", type=float, default=0.0)
    args = parser.parse_args()

    filt_class = Filter()
    print(filt_class.U_LOW_MAG)

    # Run the simulation
    phosimDir = getPhoSimPath()

    if (args.output == ""):
        outputDir = getAoclcOutputPath()
    else:
        outputDir = args.output

    os.makedirs(outputDir, exist_ok=True)
    if (args.clobber == True):
        _eraseFolderContent(outputDir)

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir

    genFlats=True

    donutLoop = genDonuts()
    for filt_name in ['u', 'g', 'r', 'i', 'z', 'y']:

        filt_class.setFilter(FilterType.fromString(filt_name))
        magBounds = filt_class.getMagBoundary()
        starMag = magBounds[0]# Index 0 is brightest

        createCat = createPhosimCatalog()
        raShift = (args.raShift * .2) / 3600 # Convert to degrees
        decShift = (args.decShift * .2) / 3600 # Convert to degrees
        createCat.createPhosimCatalog(1, 0, [starMag], raShift, decShift,
                                      args.skyFile, numFields=1)

        donutLoop.main(phosimDir, args.numOfProc, outputDir, args.testLabel,
                       surveyFilter=filt_name, starMag=starMag,
                       isEimg=args.eimage, genFlats=genFlats, useMinDofIdx=args.minDof,
                       inputSkyFilePath=args.skyFile, m1m3ForceError=args.m1m3FErr)
