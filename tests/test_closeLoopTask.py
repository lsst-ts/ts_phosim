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
import argparse
import tempfile
import unittest

from lsst.ts.wep.Utility import CamType, FilterType
from lsst.ts.ofc.Utility import InstName

from lsst.ts.phosim.CloseLoopTask import CloseLoopTask
from lsst.ts.phosim.Utility import getModulePath


class TestCloseLoopTask(unittest.TestCase):
    """Test the CloseLoopTask class."""

    def setUp(self):

        self.closeLoopTask = CloseLoopTask()

        rootTestDir = os.path.join(getModulePath(), "tests")
        self.testDir = tempfile.TemporaryDirectory(dir=rootTestDir)

    def tearDown(self):

        self.testDir.cleanup()

    def testConfigSkySimWithError(self):

        self.assertRaises(ValueError, self.closeLoopTask.configSkySim, "NoThisInstName")

    def testConfigSkySimNoSkyFileLSSTFAM(self):

        self.closeLoopTask.configSkySim(InstName.LSSTFAM)

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 189)

    def testConfigSkySimNoSkyFileCOMCAM(self):

        self.closeLoopTask.configSkySim(InstName.COMCAM)

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 9)

    def testConfigSkySimNoSkyFileLSST(self):

        self.closeLoopTask.configSkySim(InstName.LSST)

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 4)

    def testConfigSkySimWithSkyFile(self):

        testSkyFile = os.path.join(
            getModulePath(), "tests", "testData", "sky", "skyComCam.txt"
        )
        self.closeLoopTask.configSkySim(InstName.COMCAM, pathSkyFile=testSkyFile)

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 9)
        self.assertEqual(skySim.getStarMag()[0], 15)

    def testConfigWepCalc(self):

        filterType = FilterType.R
        pathIsrDir = self._makeIsrDir()
        boresight = [1.2, 2.3]
        rotCamInDeg = 30
        self.closeLoopTask.configWepCalc(
            CamType.LsstFamCam,
            pathIsrDir,
            filterType,
            boresight,
            rotCamInDeg,
            useEimg=True,
        )

        wepCalc = self.closeLoopTask.getWepCalc()
        self.assertEqual(wepCalc.getIsrDir(), pathIsrDir)
        self.assertEqual(wepCalc.getFilter(), filterType)
        self.assertEqual(wepCalc.getBoresight(), tuple(boresight))
        self.assertEqual(wepCalc.getRotAng(), rotCamInDeg)

        setting = wepCalc.getSettingFile()
        self.assertEqual(setting.getSetting("imageType"), "eimage")

    def _makeIsrDir(self):

        isrDirName = "inputTest"
        isrDir = self.closeLoopTask.createIsrDir(
            self.testDir.name, isrDirName=isrDirName
        )

        return isrDir

    def testConfigOfcCalc(self):

        instName = InstName.COMCAM
        filterType = FilterType.R
        rotAngInDeg = 30
        self.closeLoopTask.configOfcCalc(instName, filterType, rotAngInDeg)

        ofcCalc = self.closeLoopTask.getOfcCalc()
        self.assertEqual(ofcCalc.getFilter(), filterType)
        self.assertEqual(ofcCalc.getRotAng(), rotAngInDeg)

    def testConfigPhosimCmpt(self):

        # Set the environment variable of phosim path
        PHOSIMPATH = "/path/to/phosim"
        os.environ["PHOSIMPATH"] = PHOSIMPATH

        # Configure the PhoSim component
        filterType = FilterType.R
        rotAngInDeg = 30
        m1m3ForceError = 0.07
        numPro = 2
        boresight = [1.1, 2.3]
        zAngleInDeg = 31.2
        seedNum = 7
        self.closeLoopTask.configPhosimCmpt(
            filterType,
            rotAngInDeg,
            m1m3ForceError,
            numPro,
            boresight=boresight,
            zAngleInDeg=zAngleInDeg,
            seedNum=seedNum,
        )
        phosimCmpt = self.closeLoopTask.getPhosimCmpt()

        self.assertEqual(phosimCmpt.getSeedNum(), seedNum)

        setting = phosimCmpt.getSettingFile()
        self.assertEqual(setting.getSetting("numPro"), 2)

        tele = phosimCmpt.getTele()
        surveyParam = tele.getSurveyParam()
        self.assertEqual(surveyParam["filterType"], filterType)
        self.assertEqual(surveyParam["rotAngInDeg"], rotAngInDeg)
        self.assertEqual(surveyParam["boresight"], tuple(boresight))
        self.assertEqual(surveyParam["zAngleInDeg"], zAngleInDeg)

        # Pop out the environment variable of phosim path
        os.environ.pop("PHOSIMPATH")

    def testGetCamTypeAndInstNameComCam(self):

        camType, instName = self.closeLoopTask.getCamTypeAndInstName("comcam")
        self.assertEqual(camType, CamType.ComCam)
        self.assertEqual(instName, InstName.COMCAM)

    def testGetCamTypeAndInstNameLsstFam(self):

        camType, instName = self.closeLoopTask.getCamTypeAndInstName("lsstfam")
        self.assertEqual(camType, CamType.LsstFamCam)
        self.assertEqual(instName, InstName.LSSTFAM)

    def testGetCamTypeAndInstNameErr(self):

        self.assertRaises(
            ValueError, self.closeLoopTask.getCamTypeAndInstName, "noThisInst"
        )

    def testGetFilterTypeRef(self):

        filterType = self.closeLoopTask.getFilterType("ref")
        self.assertEqual(filterType, FilterType.REF)

    def testGetFilterTypeU(self):

        filterType = self.closeLoopTask.getFilterType("u")
        self.assertEqual(filterType, FilterType.U)

    def testGetFilterTypeG(self):

        filterType = self.closeLoopTask.getFilterType("g")
        self.assertEqual(filterType, FilterType.G)

    def testGetFilterTypeR(self):

        filterType = self.closeLoopTask.getFilterType("r")
        self.assertEqual(filterType, FilterType.R)

    def testGetFilterTypeI(self):

        filterType = self.closeLoopTask.getFilterType("i")
        self.assertEqual(filterType, FilterType.I)

    def testGetFilterTypeZ(self):

        filterType = self.closeLoopTask.getFilterType("z")
        self.assertEqual(filterType, FilterType.Z)

    def testGetFilterTypeY(self):

        filterType = self.closeLoopTask.getFilterType("y")
        self.assertEqual(filterType, FilterType.Y)

    def testGetFilterTypeErr(self):

        filterType = self.closeLoopTask.getFilterType("y")
        self.assertEqual(filterType, FilterType.Y)

    def testCheckAndCreateBaseOutputDir(self):

        self.assertRaises(
            ValueError, self.closeLoopTask.getFilterType, "noThisFilterType"
        )

    def testSetDefaultParser(self):

        parser = argparse.ArgumentParser()
        parser = CloseLoopTask.setDefaultParser(parser)

        args = parser.parse_known_args()[0]
        self.assertEqual(args.inst, "comcam")
        self.assertEqual(args.filterType, "ref")
        self.assertEqual(args.rotCam, 0.0)
        self.assertEqual(args.m1m3FErr, 0.05)
        self.assertEqual(args.numOfProc, 1)
        self.assertEqual(args.iterNum, 5)
        self.assertEqual(args.output, "")
        self.assertEqual(args.clobber, False)

    def testSetImgParser(self):

        parser = argparse.ArgumentParser()
        parser = CloseLoopTask.setImgParser(parser)

        argsToTest = ["--boresightDeg", "1.2", "2.3"]

        args = parser.parse_known_args(args=argsToTest)[0]
        self.assertEqual(args.boresightDeg, [1.2, 2.3])
        self.assertEqual(args.skyFile, "")
        self.assertEqual(args.eimage, False)

    def testGetSensorNameListOfFieldsComCam(self):

        sensorNameList = self.closeLoopTask.getSensorNameListOfFields(InstName.COMCAM)
        self.assertEqual(len(sensorNameList), 9)

        sensorNameListAns = [
            "R22_S00",
            "R22_S01",
            "R22_S02",
            "R22_S10",
            "R22_S11",
            "R22_S12",
            "R22_S20",
            "R22_S21",
            "R22_S22",
        ]
        self.assertEqual(sensorNameList, sensorNameListAns)

    def testGetSensorNameListOfFieldsLsstFam(self):

        sensorNameList = self.closeLoopTask.getSensorNameListOfFields(InstName.LSSTFAM)
        self.assertEqual(len(sensorNameList), 189)
        self.assertEqual(sensorNameList[0], "R01_S00")
        self.assertEqual(sensorNameList[3], "R02_S00")
        self.assertEqual(sensorNameList[5], "R02_S02")
        self.assertEqual(sensorNameList[9], "R01_S10")
        self.assertEqual(sensorNameList[13], "R02_S11")
        self.assertEqual(sensorNameList[188], "R43_S22")

    def testEraseDirectoryContent(self):

        # Make the temporary directory
        tempDir = os.path.join(self.testDir.name, "tempDir")
        os.mkdir(tempDir)
        files = os.listdir(self.testDir.name)
        self.assertEqual(len(files), 1)

        # Try to erase the content
        self.closeLoopTask.eraseDirectoryContent(self.testDir.name)

        files = os.listdir(self.testDir.name)
        self.assertEqual(len(files), 0)

    def testCheckBoresight(self):

        self.assertRaises(ValueError, self.closeLoopTask.checkBoresight, [1, 2], "")
        self.assertRaises(ValueError, self.closeLoopTask.checkBoresight, [-1, 0], "")
        self.assertRaises(ValueError, self.closeLoopTask.checkBoresight, [361, 0], "")
        self.assertRaises(ValueError, self.closeLoopTask.checkBoresight, [0, -91], "")
        self.assertRaises(ValueError, self.closeLoopTask.checkBoresight, [0, 91], "")

    def testCreateIsrDir(self):

        isrDir = self._makeIsrDir()
        self.assertTrue(os.path.exists(isrDir))

    def testMakeCalibsLSST(self):

        fakeFlatDir = self.closeLoopTask.makeCalibs(InstName.LSST, self.testDir.name)

        self.assertTrue(os.path.exists(fakeFlatDir))

        files = os.listdir(fakeFlatDir)
        self.assertEqual(len(files), 24)

    def testMakeCalibsComCam(self):

        fakeFlatDir = self.closeLoopTask.makeCalibs(InstName.COMCAM, self.testDir.name)

        self.assertTrue(os.path.exists(fakeFlatDir))

        files = os.listdir(fakeFlatDir)
        self.assertEqual(len(files), 54)

    def testAssignImgType(self):

        self.assertFalse(self.closeLoopTask.useCcdImg())

        self.closeLoopTask.assignImgType(False)
        self.assertTrue(self.closeLoopTask.useCcdImg())

        self.closeLoopTask.assignImgType(True)
        self.assertTrue(self.closeLoopTask.useCcdImg())

        self.closeLoopTask.assignImgType(None)
        self.assertFalse(self.closeLoopTask.useCcdImg())

    def testUseCcdImg(self):

        self.assertFalse(self.closeLoopTask.useCcdImg())

        self.closeLoopTask.assignImgType(False)
        self.assertTrue(self.closeLoopTask.useCcdImg())

    def testSetWepCalcWithSkyInfo(self):

        self.closeLoopTask.configSkySim(InstName.LSST)

        pathIsrDir = self._makeIsrDir()
        self.closeLoopTask.configWepCalc(
            CamType.LsstFamCam, pathIsrDir, FilterType.R, [0, 0], 0,
        )

        outputSkyInfoFilePath = self.closeLoopTask.setWepCalcWithSkyInfo(
            self.testDir.name
        )
        self.assertTrue(os.path.exists(outputSkyInfoFilePath))

        wepCalc = self.closeLoopTask.getWepCalc()
        self.assertEqual(wepCalc.getSkyFile(), outputSkyInfoFilePath)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
