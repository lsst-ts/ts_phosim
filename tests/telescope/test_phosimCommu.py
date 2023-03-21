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

from lsst.ts.wep.utility import FilterType

from lsst.ts.phosim.telescope.PhosimCommu import PhosimCommu

from lsst.ts.phosim.utils.Utility import getModulePath, SurfaceType


class TestPhosimCommu(unittest.TestCase):
    """Test the PhosimCommu class."""

    def setUp(self):
        self.phosimCom = PhosimCommu()

    def testSetAndGetPhoSimDir(self):
        phosimDir = "NotPhoSimDir"
        self.phosimCom.setPhoSimDir(phosimDir)
        self.assertEqual(self.phosimCom.getPhoSimDir(), phosimDir)

    def testGetFilterId(self):
        phosimFilterID = self.phosimCom.getFilterId(FilterType.R)
        self.assertEqual(phosimFilterID, 2)

    def testGetFilterIdWithFilterRef(self):
        self.assertRaises(ValueError, self.phosimCom.getFilterId, FilterType.REF)

    def testGetSurfaceId(self):
        surfaceID = self.phosimCom.getSurfaceId(SurfaceType.M3)
        self.assertEqual(surfaceID, 2)

    def testDoDofPert(self):
        dofInUm = np.zeros(50)
        content = self.phosimCom.doDofPert(dofInUm)
        self.assertEqual(len(content.split("\n")), self.phosimCom.DOF_NUM + 1)
        self.assertRaises(ValueError, self.phosimCom.doDofPert, np.zeros(49))

    def testDoSurfPert(self):
        surfId = 1
        zkInMm = [1]
        content = self.phosimCom.doSurfPert(surfId, zkInMm)
        ansContent = "izernike 1 0 1 \n"
        self.assertEqual(content, ansContent)

    def testDoSurfMapPert(self):
        surfId = 1
        surfFilePath = "temp.txt"
        surfFilePath = os.path.abspath(surfFilePath)
        relScale = 1
        content = self.phosimCom.doSurfMapPert(surfId, surfFilePath, relScale)
        ansContent = "surfacemap 1 %s 1 \n" % surfFilePath
        self.assertEqual(content, ansContent)

    def testDoCameraConfig(self):
        content = self.phosimCom.doCameraConfig(guidSensorOn=True)
        ansContent = "camconfig 4 \n"
        self.assertEqual(content, ansContent)

    def testDoSurfLink(self):
        linkSurfId1 = 1
        linkSurfId2 = 2
        content = self.phosimCom.doSurfLink(linkSurfId1, linkSurfId2)
        ansContent = "surfacelink 1 2 \n"
        self.assertEqual(content, ansContent)

    def testGenerateOpd(self):
        opdId = 0
        fieldXInDeg = 1.0
        fieldYInDeg = 2.0
        wavelengthInNm = 500
        content = self.phosimCom.generateOpd(
            opdId, fieldXInDeg, fieldYInDeg, wavelengthInNm
        )
        ansContent = "opd  0\t 1.000000\t 2.000000 500.0 \n"
        self.assertEqual(content, ansContent)

    def testGenerateStar(self):
        starId = 0
        ra = 1.0
        dec = 1.0
        magNorm = 2.0
        sedName = "flat.txt"
        content = self.phosimCom.generateStar(starId, ra, dec, magNorm, sedName)
        ansContent = "object  0\t 1.000000\t 1.000000  2.000000 "
        ansContent += "../sky/flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 star 0.0 "
        ansContent += "none none \n"
        self.assertEqual(content, ansContent)

    def testGetStarInstance(self):
        obsId = 100
        aFilterId = 1
        content = self.phosimCom.getStarInstance(obsId, aFilterId)
        self.assertEqual(len(content.split("\n")), 8)

    def testGetOpdInstance(self):
        obsId = 100
        aFilterId = 1
        content = self.phosimCom.getOpdInstance(obsId, aFilterId)
        self.assertEqual(len(content.split("\n")), 6)

    def testWriteToFile(self):
        filePath = os.path.join(getModulePath(), "tests", "temp.inst")
        self.assertFalse(os.path.exists(filePath))

        content = "temp"
        self.phosimCom.writeToFile(filePath, content=content, mode="w")
        self.assertTrue(os.path.exists(filePath))

        with open(filePath, "r") as file:
            contentInFile = file.read()
        self.assertEqual(contentInFile, content)

        os.remove(filePath)

    def testWriteSedFile(self):
        wavelengthInNm = 300.0
        try:
            sedFilePath = self.phosimCom.writeSedFile(wavelengthInNm)
            os.remove(sedFilePath)
        except FileNotFoundError:
            print("Do not find the sed file.")

    def testGetPhoSimArgs(self):
        instFile = "temp.inst"
        sensorName = "R22_S11"
        argString = self.phosimCom.getPhoSimArgs(instFile, sensorName=sensorName)

        absFilePath = os.path.abspath(instFile)
        ansArgString = "%s -i lsst -e 1 -s %s" % (absFilePath, sensorName)

        self.assertEqual(argString, ansArgString)

    def testFunc(self):
        try:
            self.phosimCom.runPhoSim()
        except RuntimeError:
            print("Do not find PhoSim directory.")


if __name__ == "__main__":
    # Run the unit test
    unittest.main()
