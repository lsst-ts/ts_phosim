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

from lsst.ts.wep.utility import CamType, FilterType
from lsst.ts.wep.paramReader import ParamReader

from lsst.ts.phosim.CloseLoopTask import CloseLoopTask
from lsst.ts.phosim.utils.Utility import getModulePath, getAoclcOutputPath


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
        self.closeLoopTask.configSkySim("lsstfam")

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 189)

    def testConfigSkySimNoSkyFileCOMCAM(self):
        self.closeLoopTask.configSkySim("comcam")

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 9)

    def testConfigSkySimNoSkyFileLSST(self):
        self.closeLoopTask.configSkySim("lsst")

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 8)

    def testConfigSkySimWithSkyFile(self):
        testSkyFile = os.path.join(
            getModulePath(), "tests", "testData", "sky", "skyComCam.txt"
        )
        self.closeLoopTask.configSkySim("comcam", pathSkyFile=testSkyFile)

        skySim = self.closeLoopTask.getSkySim()
        self.assertEqual(len(skySim.getStarId()), 9)
        self.assertEqual(skySim.getStarMag()[0], 15)

    def _makeIsrDir(self):
        isrDirName = "inputTest"
        isrDir = self.closeLoopTask.createIsrDir(
            self.testDir.name, isrDirName=isrDirName
        )

        return isrDir

    def testConfigOfcCalc(self):
        instName = "comcam"
        self.closeLoopTask.configOfcCalc(instName)

        ofcCalc = self.closeLoopTask.getOfcCalc()
        self.assertEqual(ofcCalc.ofc_data.name, instName)

    def testConfigPhosimCmpt(self):
        # Set the environment variable of phosim path
        PHOSIMPATH = "/path/to/phosim"
        os.environ["PHOSIMPATH"] = PHOSIMPATH

        # Configure the PhoSim component
        filterType = FilterType.LSST_R
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
        self.assertEqual(instName, "comcam")

    def testGetCamTypeAndInstNameLsstFam(self):
        camType, instName = self.closeLoopTask.getCamTypeAndInstName("lsstfam")
        self.assertEqual(camType, CamType.LsstFamCam)
        self.assertEqual(instName, "lsstfam")

    def testGetCamTypeAndInstNameLsst(self):
        camType, instName = self.closeLoopTask.getCamTypeAndInstName("lsst")
        self.assertEqual(camType, CamType.LsstCam)
        self.assertEqual(instName, "lsst")

    def testGetCamTypeAndInstNameErr(self):
        self.assertRaises(
            ValueError, self.closeLoopTask.getCamTypeAndInstName, "noThisInst"
        )

    def testGetFilterTypeRef(self):
        filterType = self.closeLoopTask.getFilterType("ref")
        self.assertEqual(filterType, FilterType.REF)

    def testGetFilterTypeU(self):
        filterType = self.closeLoopTask.getFilterType("u")
        self.assertEqual(filterType, FilterType.LSST_U)

    def testGetFilterTypeG(self):
        filterType = self.closeLoopTask.getFilterType("g")
        self.assertEqual(filterType, FilterType.LSST_G)

    def testGetFilterTypeR(self):
        filterType = self.closeLoopTask.getFilterType("r")
        self.assertEqual(filterType, FilterType.LSST_R)

    def testGetFilterTypeI(self):
        filterType = self.closeLoopTask.getFilterType("i")
        self.assertEqual(filterType, FilterType.LSST_I)

    def testGetFilterTypeZ(self):
        filterType = self.closeLoopTask.getFilterType("z")
        self.assertEqual(filterType, FilterType.LSST_Z)

    def testGetFilterTypeY(self):
        filterType = self.closeLoopTask.getFilterType("y")
        self.assertEqual(filterType, FilterType.LSST_Y)

    def testGetFilterTypeErr(self):
        self.assertRaises(
            ValueError, self.closeLoopTask.getFilterType, "noThisFilterType"
        )

    def testGetMagLimits(self):
        # check that each magLimits dictionary contains both low
        # and high values
        for filterTypeName in "ugrizy":
            magLimits = self.closeLoopTask.getMagLimits(filterTypeName)
            self.assertCountEqual(magLimits.keys(), ["low", "high"])

        # check that incorrect filter name throws an exception
        self.assertRaises(ValueError, self.closeLoopTask.getMagLimits, "x")

    def testCheckAndCreateBaseOutputDir(self):
        # first, check that the output dir is created in the AOCLC output
        # dir if no output dir is provided
        self.closeLoopTask.checkAndCreateBaseOutputDir("")
        self.assertTrue(os.path.exists(getAoclcOutputPath()))

        # second, check that the output dir is created
        # where the name is given
        baseOutputDir = os.path.join(self.testDir.name, "testBaseOutputDir")
        self.assertFalse(os.path.exists(baseOutputDir))
        self.closeLoopTask.checkAndCreateBaseOutputDir(baseOutputDir)
        self.assertTrue(os.path.exists(baseOutputDir))

    def testSetDefaultParser(self):
        parser = argparse.ArgumentParser()
        parser = CloseLoopTask.setDefaultParser(parser)

        args = parser.parse_known_args()[0]
        self.assertEqual(args.inst, "comcam")
        self.assertEqual(args.filterType, "")
        self.assertEqual(args.rotCam, 0.0)
        self.assertEqual(args.m1m3FErr, 0.05)
        self.assertEqual(args.numOfProc, 1)
        self.assertEqual(args.iterNum, 5)
        self.assertEqual(args.output, "")
        self.assertEqual(args.clobber, False)
        self.assertEqual(args.pipelineFile, "")

    def testSetImgParser(self):
        parser = argparse.ArgumentParser()
        parser = CloseLoopTask.setImgParser(parser)

        argsToTest = ["--boresightDeg", "1.2", "2.3"]

        args = parser.parse_known_args(args=argsToTest)[0]
        self.assertEqual(args.boresightDeg, [1.2, 2.3])
        self.assertEqual(args.skyFile, "")
        self.assertEqual(args.eimage, False)

    def testGetSensorNameListOfFieldsComCam(self):
        sensorNameList = self.closeLoopTask.getSensorNameListOfFields("comcam")
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
        sensorNameList = self.closeLoopTask.getSensorNameListOfFields("lsstfam")
        self.assertEqual(len(sensorNameList), 189)
        self.assertEqual(sensorNameList[0], "R01_S00")
        self.assertEqual(sensorNameList[3], "R01_S10")
        self.assertEqual(sensorNameList[5], "R01_S12")
        self.assertEqual(sensorNameList[9], "R02_S00")
        self.assertEqual(sensorNameList[13], "R02_S11")
        self.assertEqual(sensorNameList[188], "R43_S22")

        # Test the wavefront detector is not in the list
        self.assertTrue("R00_SW0" not in sensorNameList)
        # Test guider detector not in the list
        self.assertTrue("R00_SG0" not in sensorNameList)

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
        fakeFlatDir = self.closeLoopTask.makeCalibs("lsst", self.testDir.name)

        self.assertTrue(os.path.exists(fakeFlatDir))

        files = os.listdir(fakeFlatDir)
        self.assertEqual(len(files), 48)

    def testMakeCalibsComCam(self):
        fakeFlatDir = self.closeLoopTask.makeCalibs("comcam", self.testDir.name)

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

    def testMapFilterRefToG(self):
        # test that the reference filter
        # gets mapped to g
        for filterTypeName in ["ref", ""]:
            mappedFilterName = self.closeLoopTask.mapFilterRefToG(filterTypeName)
            self.assertEqual(mappedFilterName, "g")

        # test that all other filters are
        # mapped to themselves
        for filterTypeName in "ugrizy":
            mappedFilterName = self.closeLoopTask.mapFilterRefToG(filterTypeName)
            self.assertEqual(mappedFilterName, filterTypeName)

    def testWriteWepConfiguration(self):
        # Check that correct instName gets written
        for instName in ["comcam", "lsst"]:
            pipelineYamlPath = os.path.join(self.testDir.name, f"{instName}_test.yaml")
            self.closeLoopTask.writeWepConfiguration(instName, pipelineYamlPath, "g")

            # test that yaml file exists
            self.assertTrue(os.path.exists(pipelineYamlPath))

            # test for correct content
            yamlFile = ParamReader(pipelineYamlPath)
            content = yamlFile.getMatContent().item()
            butlerInstName = "ComCam" if instName == "comcam" else "Cam"
            self.assertEqual(
                content["instrument"], f"lsst.obs.lsst.Lsst{butlerInstName}"
            )

        # Check that correct filter gets written
        for filterTypeName in "ugrizy":
            pipelineYamlPath = os.path.join(
                self.testDir.name, f"{instName}_{filterTypeName}_test.yaml"
            )
            self.closeLoopTask.writeWepConfiguration(
                instName, pipelineYamlPath, filterTypeName
            )

            # read the written yaml file
            yamlFile = ParamReader(pipelineYamlPath)
            content = yamlFile.getMatContent().item()

            # test that the correct content was written
            config = content["tasks"]
            self.assertTrue("isr" in config)
            self.assertTrue("generateDonutCatalogWcsTask" in config)


if __name__ == "__main__":
    # Run the unit test
    unittest.main()
