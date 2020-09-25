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


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
