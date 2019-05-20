import os
import shutil
import numpy as np
import unittest

from lsst.ts.wep.Utility import FilterType, CamType
from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData
from lsst.ts.wep.ctrlIntf.MapSensorNameAndId import MapSensorNameAndId

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.Utility import getModulePath
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt


class TestPhosimCmpt(unittest.TestCase):
    """ Test the PhosimCmpt class."""

    def setUp(self):

        self.phosimCmpt = PhosimCmpt(self.tele)

        # Set the output directories
        self.outputDir = os.path.join(getModulePath(), "tests", "tmp")
        self.outputImgDir = os.path.join(self.outputDir, "img")

        self.phosimCmpt.setOutputDir(self.outputDir)
        self.phosimCmpt.setOutputImgDir(self.outputImgDir)

        # Set the file name of analyzed OPD data
        self.zkFileName = "opd.zer"
        self.pssnFileName = "PSSN.txt"

    @classmethod
    def setUpClass(cls):
        """Only do the instantiation for one time for the slow speed."""

        cls.tele = TeleFacade()
        cls.tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
        cls.tele.setSensorOn(sciSensorOn=True, wfSensorOn=False,
                             guidSensorOn=False)
        cls.tele.setInstName(CamType.ComCam)

        # Set the survey parameters
        obsId = 9006000
        filterType = FilterType.REF
        boresight = (0.2, 0.3)
        zAngleInDeg = 27.0912
        rotAngInDeg = np.rad2deg(-1.2323)
        cls.tele.setSurveyParam(obsId=obsId, filterType=filterType,
                                boresight=boresight, zAngleInDeg=zAngleInDeg,
                                rotAngInDeg=rotAngInDeg)

        # Set up the phosim directory
        phosimDir = "phosimDir"
        cls.tele.setPhoSimDir(phosimDir)

    def tearDown(self):

        self._setDefaultTeleSetting()

        shutil.rmtree(self.outputDir)

    def _setDefaultTeleSetting(self):

        self.tele.setDofInUm(np.zeros(50))

    def testGetNumOfZk(self):

        numOfZk = self.phosimCmpt.getNumOfZk()
        self.assertTrue(isinstance(numOfZk, int))
        self.assertEqual(numOfZk, 19)

    def testGetIntraFocalDirName(self):

        dirName = self.phosimCmpt.getIntraFocalDirName()
        self.assertEqual(dirName, "intra")

    def testGetExtraFocalDirName(self):

        dirName = self.phosimCmpt.getExtraFocalDirName()
        self.assertEqual(dirName, "extra")

    def testGetOutputDir(self):

        outputDir = self.phosimCmpt.getOutputDir()
        self.assertEqual(outputDir, self.outputDir)

    def testSetOutputDir(self):

        outputDir = os.path.join(self.outputDir, "testOutputDir")
        self.phosimCmpt.setOutputDir(outputDir)

        self.assertTrue(self._isDirExists(outputDir))
        self.assertEqual(self.phosimCmpt.getOutputDir(), outputDir)

    def _isDirExists(self, dirPath):

        return os.path.exists(dirPath)

    def testGetOutputImgDir(self):

        outputImgDir = self.phosimCmpt.getOutputImgDir()
        self.assertEqual(outputImgDir, self.outputImgDir)

    def testSetOutputImgDir(self):

        outputImgDir = os.path.join(self.outputDir, "testOutputImgDir")
        self.phosimCmpt.setOutputImgDir(outputImgDir)

        self.assertTrue(self._isDirExists(outputImgDir))
        self.assertEqual(self.phosimCmpt.getOutputImgDir(), outputImgDir)

    def testGetSeedNum(self):

        self.assertEqual(self.phosimCmpt.getSeedNum(), 0)

    def testSetSeedNum(self):

        seedNum = 3
        self.phosimCmpt.setSeedNum(seedNum)

        self.assertEqual(self.phosimCmpt.getSeedNum(), seedNum)

    def testGetPhosimParam(self):

        phosimParam = self.phosimCmpt.getPhosimParam()

        self.assertEqual(phosimParam["numPro"], 1)
        self.assertEqual(phosimParam["e2ADC"], 1)

    def testSetPhosimParam(self):

        numPro = 3
        e2ADC = 0
        self.phosimCmpt.setPhosimParam(numPro, e2ADC)

        phosimParam = self.phosimCmpt.getPhosimParam()

        self.assertEqual(phosimParam["numPro"], numPro)
        self.assertEqual(phosimParam["e2ADC"], e2ADC)

    def testSetPhosimParamWithWrongValue(self):

        self.assertRaises(ValueError, self.phosimCmpt.setPhosimParam, -1, 0)
        self.assertRaises(ValueError, self.phosimCmpt.setPhosimParam, 1, 2)

    def testGetComCamOpdArgsAndFilesForPhoSim(self):

        instFileName = "opd.inst"
        argString = self.phosimCmpt.getComCamOpdArgsAndFilesForPhoSim(
            instFileName=instFileName)

        self.assertTrue(isinstance(argString, str))
        self.assertEqual(self._getNumOfFileInFolder(self.outputDir), 11)

        instFilePath = os.path.join(self.outputDir, instFileName)
        numOfLine = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLine, 67)

    def _getNumOfFileInFolder(self, folder):

        return len([name for name in os.listdir(folder)
                   if os.path.isfile(os.path.join(folder, name))])

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testGetStarArgsAndFilesForPhoSim(self):

        skySim = self._addSglStarToSkySim()

        instFileName = "star.inst"
        argString = self.phosimCmpt.getStarArgsAndFilesForPhoSim(
            skySim, instFileName=instFileName)

        self.assertTrue(isinstance(argString, str))
        self.assertEqual(self._getNumOfFileInFolder(self.outputDir), 11)

        instFilePath = os.path.join(self.outputDir, instFileName)
        numOfLine = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLine, 63)

    def _addSglStarToSkySim(self):

        skySim = SkySim()
        skySim.addStarByRaDecInDeg(0, 0.1, 0.2, 5.0)

        return skySim

    def testSaveDofInUmFileForNextIter(self):

        dofInUm = np.arange(50)
        dofInUmFileName = "dofPertInNextIter.mat"
        self.phosimCmpt.saveDofInUmFileForNextIter(
            dofInUm, dofInUmFileName=dofInUmFileName)

        filePath = os.path.join(self.outputDir, dofInUmFileName)
        self.assertTrue(os.path.exists(filePath))

        data = np.loadtxt(filePath)
        delta = np.sum(np.abs(dofInUm - data))
        self.assertEqual(delta, 0)

    def testAnalyzeComCamOpdData(self):

        self._analyzeComCamOpdData()

        zkFilePath = os.path.join(self.outputImgDir, self.zkFileName)
        pssnFilePath = os.path.join(self.outputImgDir, self.pssnFileName)
        self.assertTrue(os.path.exists(zkFilePath))
        self.assertTrue(os.path.exists(pssnFilePath))

        zk = np.loadtxt(zkFilePath)
        ansZkFilePath = os.path.join(self._getOpdFileDirOfComCam(),
                                     "sim7_iter0_opd.zer")
        ansZk = np.loadtxt(ansZkFilePath)

        delta = np.sum(np.abs(zk - ansZk[:, 3:]))
        self.assertLess(delta, 1e-10)

        pssnData = np.loadtxt(pssnFilePath)
        pssn = pssnData[0, :]
        ansPssnFilePath = os.path.join(self._getOpdFileDirOfComCam(),
                                       "sim7_iter0_PSSN.txt")
        ansPssnData = np.loadtxt(ansPssnFilePath)
        ansPssn = ansPssnData[0, :]

        delta = np.sum(np.abs(pssn - ansPssn))
        self.assertLess(delta, 1e-10)

    def _analyzeComCamOpdData(self):

        self._copyOpdToImgDirFromTestData()
        self.phosimCmpt.analyzeComCamOpdData(
            zkFileName=self.zkFileName, pssnFileName=self.pssnFileName)

    def _copyOpdToImgDirFromTestData(self):

        opdFileDir = self._getOpdFileDirOfComCam()
        opdFileList = self.phosimCmpt._getOpdFileInDir(opdFileDir)
        for opdFile in opdFileList:
            shutil.copy2(opdFile, self.outputImgDir)

    def _getOpdFileDirOfComCam(self):

        opdFileDir = os.path.join(getModulePath(), "tests", "testData",
                                  "comcamOpdFile", "iter0")

        return opdFileDir

    def testMapOpdDataToListOfWfErr(self):

        self._analyzeComCamOpdData()

        refSensorNameList = ["R22_S00", "R22_S01", "R22_S02", "R22_S10",
                             "R22_S11", "R22_S12", "R22_S20", "R22_S21",
                             "R22_S22"]
        listOfWfErr = self.phosimCmpt.mapOpdDataToListOfWfErr(
            self.zkFileName, refSensorNameList)

        self.assertEqual(len(listOfWfErr), len(refSensorNameList))

        opdZk = self.phosimCmpt._getZkFromFile(self.zkFileName)
        mapSensorNameAndId = MapSensorNameAndId()
        for wfErr, refSensorName, zk in zip(listOfWfErr, refSensorNameList, opdZk):

            sensorId = wfErr.getSensorId()
            sensorNameList = mapSensorNameAndId.mapSensorIdToName(sensorId)[0]
            self.assertEqual(sensorNameList[0], refSensorName)

            zkInWfErr = wfErr.getAnnularZernikePoly()
            delta = np.sum(np.abs(zkInWfErr - zk))
            self.assertEqual(delta, 0)

    def testGetZkFromFile(self):

        self._analyzeComCamOpdData()

        # The correctness of values have been tested at the test case of
        # testAnalyzeComCamOpdData
        zk = self.phosimCmpt._getZkFromFile(self.zkFileName)
        self.assertEqual(zk.shape, (9, 19))

    def testGetOpdPssnFromFile(self):

        self._analyzeComCamOpdData()

        # The correctness of values have been tested at the test case of
        # testAnalyzeComCamOpdData
        pssn = self.phosimCmpt.getOpdPssnFromFile(self.pssnFileName)
        self.assertEqual(len(pssn), 9)

    def testGetOpdGqEffFwhmFromFile(self):

        self._analyzeComCamOpdData()

        gqEffFwhm = self.phosimCmpt.getOpdGqEffFwhmFromFile(
            self.pssnFileName)
        self.assertAlmostEqual(gqEffFwhm, 0.5534, places=3)

    def testGetOpdMetr(self):

        metr = self.phosimCmpt.getOpdMetr()
        self.assertTrue(isinstance(metr, OpdMetrology))

    def testAddOpdFieldXYbyDeg(self):

        fieldXInDegree = 1.2
        fieldYInDegree = 1.3
        self.phosimCmpt.addOpdFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

        metr = self.phosimCmpt.getOpdMetr()
        self.assertEqual(len(metr.fieldX), 1)
        self.assertEqual(len(metr.fieldY), 1)

        self.assertEqual(metr.fieldX[0], fieldXInDegree)
        self.assertEqual(metr.fieldY[0], fieldYInDegree)

    def testGetDofInUm(self):

        dofInUm = self.phosimCmpt.getDofInUm()
        self.assertEqual(len(dofInUm), 50)
        self.assertEqual(np.sum(np.abs(dofInUm)), 0)

    def testAccDofInUm(self):

        accDofInUm = np.arange(50)
        repeatTimes = 2
        for ii in range(repeatTimes):
            self.phosimCmpt.accDofInUm(accDofInUm)

        dofInUm = self.phosimCmpt.getDofInUm()
        delta = np.sum(np.abs(dofInUm - repeatTimes * accDofInUm))
        self.assertEqual(delta, 0)

    def testSetDofInUm(self):

        settingDofInUm = np.arange(50)
        self.phosimCmpt.setDofInUm(settingDofInUm)

        dofInUm = self.phosimCmpt.getDofInUm()
        delta = np.sum(np.abs(dofInUm - settingDofInUm))
        self.assertEqual(delta, 0)

    def testGetComCamStarArgsAndFilesForPhoSim(self):

        extraObsId = 9005000
        intraObsId = 9005001
        skySim = self._addSglStarToSkySim()
        argStringList = self.phosimCmpt.getComCamStarArgsAndFilesForPhoSim(
            extraObsId, intraObsId, skySim, simSeed=1000,
            cmdSettingFileName="starDefault.cmd",
            instSettingFileName="starSingleExp.inst")

        self.assertEqual(len(argStringList), 2)
        self.assertTrue(isinstance(argStringList[0], str))
        self.assertEqual(self._getNumOfFileInFolder(self.outputDir), 12)

        instFilePath = os.path.join(self.phosimCmpt.getOutputDir(),
                                    "starExtra.inst")
        numOfLine = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLine, 63)

    def testRepackageComCamImgFromPhoSim(self):

        self._copyComCamAmpFiles()

        intraFileFolderPath = os.path.join(
            self.phosimCmpt.getOutputImgDir(),
            self.phosimCmpt.getIntraFocalDirName())
        intraFileNum = self._getNumOfFileInFolder(intraFileFolderPath)
        self.assertEqual(intraFileNum, 16)

        self.phosimCmpt.repackageComCamImgFromPhoSim()

        intraFileNum = self._getNumOfFileInFolder(intraFileFolderPath)
        self.assertEqual(intraFileNum, 1)

        extraFileFolderPath = os.path.join(
            self.phosimCmpt.getOutputImgDir(),
            self.phosimCmpt.getExtraFocalDirName())
        extraFileNum = self._getNumOfFileInFolder(extraFileFolderPath)
        self.assertEqual(extraFileNum, 1)

    def _copyComCamAmpFiles(self):

        intraFocalDirName = self.phosimCmpt.getIntraFocalDirName()
        extraFocalDirName = self.phosimCmpt.getExtraFocalDirName()
        for imgType in (intraFocalDirName, extraFocalDirName):
            ampDirPath = os.path.join(getModulePath(), "tests", "testData",
                                      "comcamAmpPhosimData", imgType)
            dst = os.path.join(self.phosimCmpt.getOutputImgDir(), imgType)
            shutil.copytree(ampDirPath, dst)

    def testReorderAndSaveWfErrFile(self):

        listOfWfErr = self._prapareListOfWfErr()

        refSensorNameList = ["R00_S12", "R00_S21", "R01_S00", "R01_S01"]
        zkFileName = "testZk.zer"
        self.phosimCmpt.reorderAndSaveWfErrFile(
            listOfWfErr, refSensorNameList, zkFileName=zkFileName)

        zkFilePath = os.path.join(self.phosimCmpt.getOutputImgDir(),
                                  zkFileName)
        zkInFile = np.loadtxt(zkFilePath)

        numOfZk = self.phosimCmpt.getNumOfZk()
        self.assertEqual(zkInFile.shape,
                         (len(refSensorNameList), numOfZk))

        self.assertEqual(np.sum(zkInFile[0, :]), 0)
        self.assertEqual(np.sum(zkInFile[3, :]), 0)

        delta = np.sum(np.abs(
            zkInFile[1, :] - listOfWfErr[2].getAnnularZernikePoly()))
        self.assertLess(delta, 1e-10)

        delta = np.sum(np.abs(
            zkInFile[2, :] - listOfWfErr[1].getAnnularZernikePoly()))
        self.assertLess(delta, 1e-10)

    def _prapareListOfWfErr(self):

        numOfZk = self.phosimCmpt.getNumOfZk()

        sensorIdList = [2, 3, 1]
        listOfWfErr = []
        for sensorId in sensorIdList:
            sensorWavefrontData = SensorWavefrontData()
            sensorWavefrontData.setSensorId(sensorId)

            wfErr = np.random.rand(numOfZk)
            sensorWavefrontData.setAnnularZernikePoly(wfErr)

            listOfWfErr.append(sensorWavefrontData)

        return listOfWfErr

    def testGetWfErrValuesAndStackToMatrix(self):

        wfErrMap, wfsValueMatrix = self._prepareWfErrMap()
        valueMatrix = \
            self.phosimCmpt._getWfErrValuesAndStackToMatrix(wfErrMap)

        delta = np.sum(np.abs(valueMatrix - wfsValueMatrix))
        self.assertEqual(delta, 0)

    def _prepareWfErrMap(self):

        sensorNameList = ["c", "b", "a"]
        numOfZk = self.phosimCmpt.getNumOfZk()
        wfsValueMatrix = np.random.rand(len(sensorNameList), numOfZk)
        wfErrMap = dict()
        for idx, sensorName in enumerate(sensorNameList):
            wfErrMap[sensorName] = wfsValueMatrix[idx, :]

        return wfErrMap, wfsValueMatrix


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
