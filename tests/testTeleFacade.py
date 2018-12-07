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

        configFilePath = os.path.join(getModulePath(), "configData",
                                      "telescopeConfig", "GT.inst")
        self.tele = TeleFacade(configFilePath=configFilePath)


        # # Set the subsystem data directory
        # self.camDataDir = os.path.join(getModulePath(), "configData", "camera")
        # self.M1M3dataDir = os.path.join(getModulePath(), "configData", "M1M3")
        # self.M2dataDir = os.path.join(getModulePath(), "configData", "M2")

        # # Set the output dir
        # self.outputDir = os.path.join(getModulePath(), "output", "temp")
        # os.makedirs(self.outputDir)

        # # Set the command setting file
        # self.starCmdSettingFile = os.path.join(getModulePath(), "configData",
        #                                        "cmdFile", "starDefault.cmd")

        # # Set the instance setting file
        # self.starInstSettingFile = os.path.join(getModulePath(), "configData",
        #                                         "instFile", "starDefault.inst")
        # self.opdInstSettingFile = os.path.join(getModulePath(), "configData",
        #                                        "instFile", "opdDefault.inst")

        # # Instantiate the needed objects
        # self.metr = OpdMetrology()
        # self.skySim = SkySim()

    def testSetConfigFile(self):

        configFilePath = "NotConfigFilePath"
        self.tele.setConfigFile(configFilePath)
        self.assertEqual(self.tele.configFile, configFilePath)

    def testGetConfigValue(self):

        varName = "M1M3TxGrad"
        value = self.tele.getConfigValue(varName)
        self.assertEqual(value, -0.0894)



    # def testFunc(self):
        
    #     # Set the configuration file path
    #     tele.setConfigFile(self.configFilePath)
    #     self.assertEqual(tele.configFile, self.configFilePath)

    #     instFilePath = "temp.inst"
    #     argString = tele.getPhoSimArgs(instFilePath)
    #     ansArgString = "%s -i lsst -e 1" % os.path.abspath(instFilePath)
    #     self.assertEqual(argString, ansArgString)

    #     instruFile = "comcam10"
    #     tele.setInstName(instruFile)
    #     self.assertEqual(tele.instName, "comcam")
    #     self.assertEqual(tele.defocalDisInMm, 1.0)

    #     dofInUm = np.random.rand(50)
    #     tele.setDofInUm(dofInUm)
    #     self.assertEqual(np.sum(np.abs(tele.dofInUm-dofInUm)), 0)

    #     tele.setDofInUm(np.zeros(50))
    #     tele.accDofInUm(dofInUm)
    #     self.assertEqual(np.sum(np.abs(tele.dofInUm-dofInUm)), 0)

    #     tele.setSubSysConfigFile(camDataDir=self.camDataDir,
    #                              M1M3dataDir=self.M1M3dataDir, 
    #                              M2dataDir=self.M2dataDir)
    #     self.assertEqual(tele.cam.camDataDir, self.camDataDir)
    #     self.assertEqual(tele.phoSimCommu.phosimDir, "")

    #     dofFilePath = tele.writeAccDofFile(self.outputDir)
    #     self.assertLess(np.sum(np.abs(np.loadtxt(dofFilePath)-dofInUm)), 1e-7)
    #     os.remove(dofFilePath)

    #     zAngleInDeg = 27.0912
    #     rotAngInDeg = -1.2323/np.pi*180.0
    #     iSim = 6
    #     pertCmdFilePath = tele.writePertBaseOnConfigFile(self.outputDir, zAngleInDeg=zAngleInDeg, 
    #                                        rotAngInDeg=rotAngInDeg, seedNum=iSim, saveResMapFig=False)
    #     cmdFile = open(pertCmdFilePath, "r")
    #     lines = cmdFile.readlines()
    #     cmdFile.close()
    #     self.assertEqual(len(lines), 256)

    #     cmdFilePath = tele.writeCmdFile(self.outputDir, cmdSettingFile=self.starCmdSettingFile, 
    #                                     pertFilePath=pertCmdFilePath, cmdFileName="star.cmd")

    #     starCmdFile = open(cmdFilePath, "r")
    #     lines = starCmdFile.readlines()
    #     starCmdFile.close()
    #     self.assertEqual(len(lines), 267)

    #     metr.addFieldXYbyDeg(0, 0)
    #     obsId = 9006000
    #     # aFilter = "g"
    #     aFilter = FilterType.G
    #     wavelengthInNm = 500
    #     instFilePath = tele.writeOpdInstFile(self.outputDir, metr, obsId, aFilter, wavelengthInNm, 
    #                                             instSettingFile=self.opdInstSettingFile)
    #     opdInstFile = open(instFilePath, "r")
    #     lines = opdInstFile.readlines()
    #     opdInstFile.close()
    #     self.assertEqual(len(lines), 55)

    #     skySim.addStarByRaDecInDeg(0, 1.0, 1.0, 17.0)
    #     boresight = (0.2, 0.3)
    #     instFilePath = tele.writeStarInstFile(self.outputDir, skySim, obsId, aFilter, boresight, 
    #                                             wfSensorOn=True, instSettingFile=self.starInstSettingFile)

    #     starInstFile = open(instFilePath, "r")
    #     lines = starInstFile.readlines()
    #     starInstFile.close()
    #     self.assertEqual(len(lines), 62)

    #     shutil.rmtree(self.outputDir)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
