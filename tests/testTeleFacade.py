import os
import shutil
import numpy as np
import unittest

from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.TeleFacade import TeleFacade

from lsst.ts.phosim.Utility import getModulePath, FilterType


class TestTeleFacade(unittest.TestCase):
    """ Test the TeleFacade class."""

    def setUp(self):

        self.configFilePath = os.path.join(getModulePath(), "configData",
                                           "telescopeConfig", "GT.inst")
        self.tele = TeleFacade(configFilePath=self.configFilePath)

        # Set the subsystem data directory
        camDataDir = os.path.join(getModulePath(), "configData", "camera")
        M1M3dataDir = os.path.join(getModulePath(), "configData", "M1M3")
        M2dataDir = os.path.join(getModulePath(), "configData", "M2")
        self.tele.setSubSysConfigDir(camDataDir=camDataDir,
                                     M1M3dataDir=M1M3dataDir,
                                     M2dataDir=M2dataDir)

        # Set the output dir
        self.outputDir = os.path.join(getModulePath(), "output", "temp")
        os.makedirs(self.outputDir)

    def tearDown(self):

        shutil.rmtree(self.outputDir)

    def testSetConfigFile(self):

        configFilePath = "NotConfigFilePath"
        self.tele.setConfigFile(configFilePath)

        self.assertEqual(self.tele.configFile, configFilePath)

    def testGetConfigValue(self):

        varName = "M1M3TxGrad"
        value = self.tele.getConfigValue(varName)

        self.assertEqual(value, -0.0894)

    def testGetPhoSimArgs(self):

        instFilePath = "temp.inst"
        argString = self.tele.getPhoSimArgs(instFilePath)

        ansArgString = "%s -i lsst -e 1" % os.path.abspath(instFilePath)

        self.assertEqual(argString, ansArgString)

    def testSetInstName(self):

        instName = "comcam10"
        self.tele.setInstName(instName)

        self.assertEqual(self.tele.surveyParam["instName"], "comcam")
        self.assertEqual(self.tele.surveyParam["defocalDisInMm"], 1.0)

        instName = "temp"
        self.tele.setInstName(instName)

        self.assertEqual(self.tele.surveyParam["instName"], "temp")
        self.assertEqual(self.tele.surveyParam["defocalDisInMm"], 1.5)

        self.assertRaises(ValueError, self.tele.setInstName, "a10a")

    def testSetDofInUm(self):

        dofInUm = np.random.rand(50)
        self.tele.setDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(self.tele.dofInUm-dofInUm)), 0)

    def testAccDofInUm(self):

        dofInUm = np.random.rand(50)
        self.tele.accDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(self.tele.dofInUm-dofInUm)), 0)

    def testSetSubSysConfigDir(self):

        tele = TeleFacade(configFilePath=self.configFilePath)
        self.assertEqual(tele.cam, None)
        self.assertEqual(tele.M1M3, None)
        self.assertEqual(tele.M2, None)

        camDataDir = "NotCamDataDir"
        M1M3dataDir = "NotM1M3dataDir"
        M2dataDir = "NotM2dataDir"
        phosimDir = "NotPhosimDir"
        tele.setSubSysConfigDir(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir,
                                M2dataDir=M2dataDir, phosimDir=phosimDir)

        # Check the subsystems are instantiated after setting up the
        # configuration directory.
        self.assertNotEqual(tele.cam, None)
        self.assertNotEqual(tele.M1M3, None)
        self.assertNotEqual(tele.M2, None)

        self.assertEqual(tele.cam.camDataDir, camDataDir)
        self.assertEqual(tele.M1M3.mirrorDataDir, M1M3dataDir)
        self.assertEqual(tele.M2.mirrorDataDir, M2dataDir)
        self.assertEqual(tele.phoSimCommu.phosimDir, phosimDir)

    def testWriteAccDofFile(self):

        dofFilePath = self.tele.writeAccDofFile(self.outputDir)
        dof = np.loadtxt(dofFilePath)

        self.assertEqual(np.sum(dof), 0)
        self.assertEqual(len(dof), 50)

    def testWritePertBaseOnConfigFile(self):

        pertCmdFilePath = self._writePertBaseOnConfigFile(self.outputDir)

        numOfLineInFile = self._getNumOfLineInFile(pertCmdFilePath)
        self.assertEqual(numOfLineInFile, 256)

    def _writePertBaseOnConfigFile(self, outputDir):

        zAngleInDeg = 27.0912
        rotAngInDeg = np.rad2deg(-1.2323)
        iSim = 6

        pertCmdFilePath = self.tele.writePertBaseOnConfigFile(
                                    outputDir, zAngleInDeg=zAngleInDeg,
                                    rotAngInDeg=rotAngInDeg, seedNum=iSim,
                                    M1M3ForceError=0.05,
                                    saveResMapFig=True,
                                    pertCmdFileName="pert.cmd")

        return pertCmdFilePath

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testWriteCmdFile(self):

        starCmdSettingFile = os.path.join(getModulePath(), "configData",
                                          "cmdFile", "starDefault.cmd")

        pertCmdFilePath = self._writePertBaseOnConfigFile(self.outputDir)
        cmdFilePath = self.tele.writeCmdFile(
                                self.outputDir,
                                cmdSettingFile=starCmdSettingFile,
                                pertFilePath=pertCmdFilePath,
                                cmdFileName="star.cmd")

        numOfLineInFile = self._getNumOfLineInFile(cmdFilePath)
        self.assertEqual(numOfLineInFile, 267)

    def testGetPhoSimCamSurf(self):

        camSurfName = "L1S2zer"
        surfaceType = self.tele._getPhoSimCamSurf(camSurfName)
        self.assertEqual(surfaceType.name, "L1B")
        self.assertRaises(ValueError, self.tele._getPhoSimCamSurf,
                          "L4S2zer")

    def testWriteOpdInstFile(self):

        metr = OpdMetrology()
        metr.addFieldXYbyDeg(0, 0)
        obsId = 9006000
        aFilter = FilterType.G
        wavelengthInNm = 500
        opdInstSettingFile = os.path.join(getModulePath(), "configData",
                                          "instFile", "opdDefault.inst")

        instFilePath = self.tele.writeOpdInstFile(
                                    self.outputDir, metr, obsId, aFilter,
                                    wavelengthInNm,
                                    instSettingFile=opdInstSettingFile)

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 55)

    def testWriteStarInstFile(self):

        skySim = SkySim()
        skySim.addStarByRaDecInDeg(0, 1.0, 1.0, 17.0)

        obsId = 9006000
        aFilter = FilterType.G
        boresight = (0.2, 0.3)

        starInstSettingFile = os.path.join(getModulePath(), "configData",
                                           "instFile", "starDefault.inst")

        instFilePath = self.tele.writeStarInstFile(
                                    self.outputDir, skySim, obsId, aFilter,
                                    boresight, wfSensorOn=True,
                                    instSettingFile=starInstSettingFile)

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 62)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
