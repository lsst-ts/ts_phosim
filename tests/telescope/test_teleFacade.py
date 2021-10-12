# This file is part of ts_phosim.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import numpy as np
import unittest
import tempfile

from lsst.ts.wep.Utility import FilterType, CamType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade

from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.SkySim import SkySim

from lsst.ts.phosim.utils.Utility import getModulePath, getConfigDir


class TestTeleFacade(unittest.TestCase):
    """ Test the TeleFacade class."""

    def setUp(self):

        self.configDir = getConfigDir()
        self.modulePath = getModulePath()

        testDir = os.path.join(self.modulePath, "tests")
        self.outputDir = tempfile.TemporaryDirectory(dir=testDir)

    @classmethod
    def setUpClass(cls):
        """Only do the instantiation for one time for the slow speed."""

        cls.tele = TeleFacade()
        cls.tele.addSubSys(addCam=True, addM1M3=True, addM2=True)

        # Set the survey parameters
        obsId = 9006000
        cls.filterType = FilterType.G
        boresight = (0.2, 0.3)
        zAngleInDeg = 27.0912
        rotAngInDeg = np.rad2deg(-1.2323)
        cls.tele.setSurveyParam(
            obsId=obsId,
            filterType=cls.filterType,
            boresight=boresight,
            zAngleInDeg=zAngleInDeg,
            rotAngInDeg=rotAngInDeg,
        )

    def tearDown(self):

        self._setDefaultTeleSetting()
        self.outputDir.cleanup()

    def _setDefaultTeleSetting(self):

        self.tele.setSensorOn()
        self.tele.setDofInUm(np.zeros(50))
        self.tele.setSurveyParam(filterType=self.filterType)
        self.tele.setInstName(CamType.LsstFamCam)

    def testGetDofInUm(self):

        dofInUm = self.tele.getDofInUm()

        self.assertEqual(len(dofInUm), 50)
        self.assertEqual(np.sum(np.abs(dofInUm)), 0)

    def testGetNumOfDof(self):

        numOfDof = self.tele.getNumOfDof()
        self.assertEqual(numOfDof, 50)

    def testGetDefaultDefocalDist(self):

        defocalDist = self.tele.getDefocalDistInMm()
        self.assertEqual(defocalDist, 1.5)

    def testGetSurfGridN(self):

        surfGridN = self.tele.getSurfGridN()
        self.assertEqual(surfGridN, 200)

    def testGetCamMjd(self):

        mjd = self.tele.getCamMjd()
        self.assertEqual(mjd, 59580.0)

    def testGetRefWaveLength(self):

        refWaveLength = self.tele.getRefWaveLength()
        self.assertEqual(refWaveLength, 500)

    def testGetDefocalDisInMm(self):

        defocalDist = 1.3
        self.tele.setInstName(CamType.LsstCam, defocalDist=defocalDist)

        self.assertEqual(self.tele.getDefocalDistInMm(), defocalDist)

    def testSetSurveyParamWithCorrectInput(self):

        obsId = 100
        filterType = FilterType.U
        boresight = (10, 20)
        zAngleInDeg = 10.0
        rotAngInDeg = 11.0

        tele = TeleFacade()
        tele.setSurveyParam(
            obsId=obsId,
            filterType=filterType,
            boresight=boresight,
            zAngleInDeg=zAngleInDeg,
            rotAngInDeg=rotAngInDeg,
        )

        self.assertEqual(tele.surveyParam["obsId"], obsId)
        self.assertEqual(tele.surveyParam["filterType"], filterType)
        self.assertEqual(tele.surveyParam["boresight"], boresight)
        self.assertEqual(tele.surveyParam["zAngleInDeg"], zAngleInDeg)
        self.assertEqual(tele.surveyParam["rotAngInDeg"], rotAngInDeg)

    def testSetSurveyParamWithWrongInput(self):

        defaultObsId = self.tele.surveyParam["obsId"]

        obsId = 1.0
        self.tele.setSurveyParam(obsId=obsId)

        self.assertEqual(self.tele.surveyParam["obsId"], defaultObsId)
        self.assertEqual(self.tele.surveyParam["filterType"], FilterType.G)

    def testSetSensorOnWithCorrectInput(self):

        sciSensorOn = False
        wfSensorOn = False
        guidSensorOn = True

        self.tele.setSensorOn(
            sciSensorOn=sciSensorOn, wfSensorOn=wfSensorOn, guidSensorOn=guidSensorOn
        )
        self.assertEqual(self.tele.sensorOn["sciSensorOn"], sciSensorOn)
        self.assertEqual(self.tele.sensorOn["wfSensorOn"], wfSensorOn)
        self.assertEqual(self.tele.sensorOn["guidSensorOn"], guidSensorOn)

    def testGetPhoSimArgs(self):

        instFilePath = "temp.inst"
        argString = self.tele.getPhoSimArgs(instFilePath)

        ansArgString = "%s -i lsst -e 1" % os.path.abspath(instFilePath)

        self.assertEqual(argString, ansArgString)

    def testSetInstName(self):

        defocalDist = 1.0
        self.tele.setInstName(CamType.ComCam, defocalDist=1.0)

        self.assertEqual(self.tele.surveyParam["instName"], "comcam")
        self.assertEqual(self.tele.getDefocalDistInMm(), defocalDist)

    def testSetInstNameWithWrongCamType(self):

        self.assertRaises(ValueError, self.tele.setInstName, "wrongType")

    def testSetInstNameWithNegativeDefocalDist(self):

        self.assertRaises(ValueError, self.tele.setInstName, CamType.LsstCam, -1.5)

        self.assertRaises(ValueError, self.tele.setInstName, CamType.LsstCam, 0)

    def testSetDofInUm(self):

        dofInUm = np.random.rand(50)
        self.tele.setDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(self.tele.dofInUm - dofInUm)), 0)

    def testSetDofInUmWithWrongLeng(self):

        dofInUm = np.random.rand(45)
        self.assertRaises(ValueError, self.tele.setDofInUm, dofInUm)

    def testAccDofInUm(self):

        dofInUm = np.random.rand(50)
        self.tele.accDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(self.tele.dofInUm - dofInUm)), 0)

    def testAccDofInUmWithWrongLeng(self):

        dofInUm = np.random.rand(45)
        self.assertRaises(ValueError, self.tele.accDofInUm, dofInUm)

    def testAddSubSys(self):

        tele = TeleFacade()
        self.assertEqual(tele.cam, None)
        self.assertEqual(tele.m1m3, None)
        self.assertEqual(tele.m2, None)

        tele.addSubSys(addCam=True, addM1M3=True, addM2=True)

        self.assertNotEqual(tele.cam, None)
        self.assertNotEqual(tele.m1m3, None)
        self.assertNotEqual(tele.m2, None)

    def testWriteAccDofFile(self):

        dofFilePath = self.tele.writeAccDofFile(self.outputDir.name)
        dof = np.loadtxt(dofFilePath)

        self.assertEqual(np.sum(dof), 0)
        self.assertEqual(len(dof), 50)

    def testWritePertBaseOnConfigFile(self):

        pertCmdFilePath = self._writePertBaseOnConfigFile(self.outputDir.name)

        numOfLineInFile = self._getNumOfLineInFile(pertCmdFilePath)
        self.assertEqual(numOfLineInFile, 256)

    def _writePertBaseOnConfigFile(self, outputDir):

        iSim = 6
        pertCmdFilePath = self.tele.writePertBaseOnConfigFile(
            outputDir,
            seedNum=iSim,
            m1m3ForceError=0.05,
            saveResMapFig=True,
            pertCmdFileName="pert.cmd",
        )

        return pertCmdFilePath

    def _getNumOfLineInFile(self, filePath):

        with open(filePath, "r") as file:
            return sum(1 for line in file.readlines())

    def testWriteCmdFile(self):

        starCmdSettingFile = os.path.join(self.configDir, "cmdFile", "starDefault.cmd")

        pertCmdFilePath = self._writePertBaseOnConfigFile(self.outputDir.name)
        cmdFilePath = self.tele.writeCmdFile(
            self.outputDir.name,
            cmdSettingFile=starCmdSettingFile,
            pertFilePath=pertCmdFilePath,
            cmdFileName="star.cmd",
        )

        numOfLineInFile = self._getNumOfLineInFile(cmdFilePath)
        self.assertEqual(numOfLineInFile, 268)

    def testGetPhoSimCamSurf(self):

        camSurfName = "L1S2zer"
        surfaceType = self.tele._getPhoSimCamSurf(camSurfName)
        self.assertEqual(surfaceType.name, "L1B")
        self.assertRaises(ValueError, self.tele._getPhoSimCamSurf, "L4S2zer")

    def testWriteOpdInstFile(self):

        metr, opdInstSettingFile = self._generateOpd()

        with self.assertWarns(UserWarning):
            instFilePath = self.tele.writeOpdInstFile(
                self.outputDir.name, metr, instSettingFile=opdInstSettingFile
            )

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 59)

    def testWriteOpdInstFileWithFilterRef(self):

        self.tele.setSurveyParam(filterType=FilterType.REF)

        metr, opdInstSettingFile = self._generateOpd()
        with self.assertWarns(UserWarning):
            instFilePath = self.tele.writeOpdInstFile(
                self.outputDir.name, metr, instSettingFile=opdInstSettingFile
            )

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 59)

    def _generateOpd(self):

        metr = OpdMetrology()
        metr.addFieldXYbyDeg(0, 0)
        opdInstSettingFile = os.path.join(self.configDir, "instFile", "opdDefault.inst")

        return metr, opdInstSettingFile

    def testWriteStarInstFile(self):

        skySim, starInstSettingFile = self._generateFakeSky()

        instFilePath = self.tele.writeStarInstFile(
            self.outputDir.name, skySim, instSettingFile=starInstSettingFile
        )

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 63)

    def testWriteStarInstFileWithFilterRef(self):

        self.tele.setSurveyParam(filterType=FilterType.REF)

        skySim, starInstSettingFile = self._generateFakeSky()
        with self.assertWarns(UserWarning):
            instFilePath = self.tele.writeStarInstFile(
                self.outputDir.name, skySim, instSettingFile=starInstSettingFile
            )

        numOfLineInFile = self._getNumOfLineInFile(instFilePath)
        self.assertEqual(numOfLineInFile, 63)

    def _generateFakeSky(self):

        skySim = SkySim()
        skySim.addStarByRaDecInDeg(0, 1.0, 1.0, 17.0)

        starInstSettingFile = os.path.join(
            self.configDir, "instFile", "starDefault.inst"
        )

        return skySim, starInstSettingFile


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
