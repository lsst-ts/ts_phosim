import os
import shutil
import numpy as np
import unittest
from copy import deepcopy

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.TeleFacade import TeleFacade

from lsst.ts.phosim.Utility import getModulePath, createObservation


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

        # Set the observation parameters
        obsId = 9006000
        filterType = FilterType.G
        boresight = (0.2, 0.3)
        zAngleInDeg = 27.0912
        rotAngInDeg = np.rad2deg(-1.2323)
        mjd = 59552.3
        obs = createObservation(obsId=obsId, filterType=filterType,
                                 boresight=boresight, zAngleInDeg=zAngleInDeg,
                                 rotAngInDeg=rotAngInDeg, mjd=mjd)
        self.tele.setObservation(obs)

        # Set the output dir
        self.outputDir = os.path.join(getModulePath(), "output", "temp")
        os.makedirs(self.outputDir)

    def tearDown(self):

        shutil.rmtree(self.outputDir)

    def testSetObservation(self):
        
        obsId = 99
        obs = createObservation(obsId=obsId)
        self.tele.setObservation(obs)
        self.assertEqual(self.tele.obs.OpsimMetaData['obsHistID'], obsId)

    def testGetObservation(self):
        
        retObsId = self.tele.getObservation().OpsimMetaData['obsHistID']
        self.assertEqual(self.tele.obs.OpsimMetaData['obsHistID'], 9006000)

    def testGetDefocalDisInMm(self):

        instName = "comcam13"
        self.tele.setInstName(instName)

        self.assertEqual(self.tele.getDefocalDisInMm(), 1.3)

    def testSetSensorOnWithCorrectInput(self):

        sciSensorOn = False
        wfSensorOn = False
        guidSensorOn = True

        self.tele.setSensorOn(sciSensorOn=sciSensorOn, wfSensorOn=wfSensorOn,
                              guidSensorOn=guidSensorOn)
        self.assertEqual(self.tele.sensorOn["sciSensorOn"], sciSensorOn)
        self.assertEqual(self.tele.sensorOn["wfSensorOn"], wfSensorOn)
        self.assertEqual(self.tele.sensorOn["guidSensorOn"], guidSensorOn)

    def testSetSensorOnWithWrongInput(self):

        sciSensorOn = 0
        self.tele.setSensorOn(sciSensorOn=sciSensorOn)
        self.assertEqual(self.tele.sensorOn["sciSensorOn"], True)

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

        self.assertEqual(self.tele.instName, "comcam")
        self.assertEqual(self.tele.getDefocalDisInMm(), 1.0)

        instName = "temp"
        self.tele.setInstName(instName)

        self.assertEqual(self.tele.instName, "temp")
        self.assertEqual(self.tele.getDefocalDisInMm(), 1.5)

        self.assertRaises(ValueError, self.tele.setInstName, "a10a")

    def testSetDofInUm(self):

        dofInUm = np.random.rand(50)
        self.tele.setDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(self.tele.dofInUm - dofInUm)), 0)

    def testAccDofInUm(self):

        dofInUm = np.random.rand(50)
        self.tele.accDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(self.tele.dofInUm - dofInUm)), 0)

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

        iSim = 6
        pertCmdFilePath = self.tele.writePertBaseOnConfigFile(
            outputDir, seedNum=iSim, M1M3ForceError=0.05, saveResMapFig=True,
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
            self.outputDir, cmdSettingFile=starCmdSettingFile,
            pertFilePath=pertCmdFilePath, cmdFileName="star.cmd")

        numOfLineInFile = self._getNumOfLineInFile(cmdFilePath)
        self.assertEqual(numOfLineInFile, 267)

    def testGetPhoSimCamSurf(self):

        camSurfName = "L1S2zer"
        surfaceType = self.tele._getPhoSimCamSurf(camSurfName)
        self.assertEqual(surfaceType.name, "L1B")
        self.assertRaises(ValueError, self.tele._getPhoSimCamSurf,
                          "L4S2zer")

    def testWriteOpdInstFile(self):

        metr, opdInstSettingFile = self._generateOpd()

        instFilePath = self.tele.writeOpdInstFile(
            self.outputDir, metr, instSettingFile=opdInstSettingFile)

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 59)

    def testWriteOpdInstFileWithFilterRef(self):

        obs = deepcopy(self.tele.obs)
        obs.bandpass = FilterType.Ref.toString()
        self.tele.setObservation(obs)

        metr, opdInstSettingFile = self._generateOpd()
        instFilePath = self.tele.writeOpdInstFile(
            self.outputDir, metr, instSettingFile=opdInstSettingFile)

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 59)

    def _generateOpd(self):

        metr = OpdMetrology()
        metr.addFieldXYbyDeg(0, 0)
        opdInstOverrideFile = os.path.join(getModulePath(), "configData",
                                          "instFile", "opdDefault.inst")

        return metr, opdInstOverrideFile

    def testWriteStarInstFile(self):

        skySim, starInstOverrideFile = self._generateFakeSky()

        instFilePath = self.tele.writeStarInstFile(
            self.outputDir, skySim, instOverrideFile=starInstOverrideFile)

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 66)

    def testWriteStarInstFileWithFilterRef(self):

        obs = deepcopy(self.tele.obs)
        obs.bandpass = FilterType.Ref.toString()
        self.tele.setObservation(obs)

        skySim, starInstOverrideFile = self._generateFakeSky()
        instFilePath = self.tele.writeStarInstFile(
            self.outputDir, skySim, instOverrideFile=starInstOverrideFile)

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 66)

    def _generateFakeSky(self):

        skySim = SkySim()
        skySim.addStarByRaDecInDeg(0, 1.0, 1.0, 17.0)

        starInstOverrideFile = os.path.join(getModulePath(), "configData",
                                           "instFile", "starDefault.inst")

        return skySim, starInstOverrideFile


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
