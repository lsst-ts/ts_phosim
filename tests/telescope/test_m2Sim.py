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
import unittest
import numpy as np

from lsst.ts.phosim.telescope.M2Sim import M2Sim

from lsst.ts.phosim.utils.Utility import getModulePath


class TestM2Sim(unittest.TestCase):
    """Test the M2Sim class."""

    def setUp(self):

        self.testM2Data = os.path.join(
            getModulePath(), "tests", "testData", "testM2Func"
        )
        self.outputDir = os.path.join(getModulePath(), "output")

    @classmethod
    def setUpClass(cls):
        """Only do the instantiation for one time for the slow speed."""

        cls.m2 = M2Sim()

    def testInit(self):

        self.assertEqual(self.m2.getInnerRinM(), 0.9)
        self.assertEqual(self.m2.getOuterRinM(), 1.71)

    def testGetPrintthz(self):

        zAngleInDeg = 27.0912
        zAngleInRadian = np.deg2rad(zAngleInDeg)
        printthzInUm = self.m2.getPrintthz(zAngleInRadian)

        self.assertEqual(len(printthzInUm), 15984)
        self.assertAlmostEqual(printthzInUm[0] * 1e4, 8.65985, places=4)
        self.assertAlmostEqual(printthzInUm[1] * 1e3, 1.36859, places=4)

        printthzMax = np.max(printthzInUm)
        self.assertAlmostEqual(printthzMax * 1e2, 4.91178, places=4)
        self.assertEqual(np.where(printthzInUm == printthzMax)[0][0], 1079)

        printthzMin = np.min(printthzInUm)
        self.assertAlmostEqual(printthzMin * 1e2, -5.27582, places=4)
        self.assertEqual(np.where(printthzInUm == printthzMin)[0][0], 15731)

    def testGetTempCorr(self):

        m2TzGrad = -0.0675
        m2TrGrad = -0.1416
        tempCorrInUm = self.m2.getTempCorr(m2TzGrad, m2TrGrad)

        self.assertEqual(len(tempCorrInUm), 15984)
        self.assertAlmostEqual(tempCorrInUm[0] * 1e3, 2.94989, places=4)
        self.assertAlmostEqual(tempCorrInUm[1] * 1e3, 2.85619, places=4)

        tempCorrMax = np.max(tempCorrInUm)
        self.assertAlmostEqual(tempCorrMax * 1e3, 2.95336, places=4)
        self.assertEqual(np.where(tempCorrInUm == tempCorrMax)[0][0], 7781)

        tempCorrMin = np.min(tempCorrInUm)
        self.assertAlmostEqual(tempCorrMin * 1e3, -4.52449, places=4)
        self.assertEqual(np.where(tempCorrInUm == tempCorrMin)[0][0], 3618)

    def testGetMirrorResInMmInZemax(self):

        self._setSurfAlongZ()
        zcInMmInZemax = self.m2.getMirrorResInMmInZemax()[3]

        self.assertEqual(len(zcInMmInZemax), 28)
        self.assertAlmostEqual(zcInMmInZemax[0] * 1e6, -1.26777, places=4)
        self.assertAlmostEqual(zcInMmInZemax[1] * 1e11, 6.16316, places=4)

        zcMax = np.max(zcInMmInZemax)
        self.assertAlmostEqual(zcMax * 1e6, 3.82926, places=4)
        self.assertEqual(np.where(zcInMmInZemax == zcMax)[0][0], 3)

        zcMin = np.min(zcInMmInZemax)
        self.assertAlmostEqual(zcMin * 1e6, -1.30288, places=4)
        self.assertEqual(np.where(zcInMmInZemax == zcMin)[0][0], 20)

    def _setSurfAlongZ(self):

        zAngleInDeg = 27.0912
        zAngleInRadian = np.deg2rad(zAngleInDeg)
        printthzInUm = self.m2.getPrintthz(zAngleInRadian)

        m2TzGrad = -0.0675
        m2TrGrad = -0.1416
        tempCorrInUm = self.m2.getTempCorr(m2TzGrad, m2TrGrad)

        mirrorSurfInUm = printthzInUm + tempCorrInUm
        self.m2.setSurfAlongZ(mirrorSurfInUm)

    def testWriteMirZkAndGridResInZemax(self):

        resFile = self._writeMirZkAndGridResInZemax()
        content = np.loadtxt(resFile)

        self.assertEqual(content.shape, (41617, 4))

        deltaX = content[0, 2]
        self.assertAlmostEqual(deltaX, 17.18592965, places=7)

        self.assertAlmostEqual(content[300, 0] * 1e6, 1.50337, places=4)
        self.assertAlmostEqual(content[300, 1] * 1e9, 9.100, places=2)
        self.assertAlmostEqual(content[300, 2] * 1e9, -8.581, places=2)

        os.remove(resFile)

    def _writeMirZkAndGridResInZemax(self):

        self._setSurfAlongZ()
        resFile = os.path.join(self.outputDir, "M2res.txt")
        self.m2.writeMirZkAndGridResInZemax(resFile=resFile)

        return resFile

    def testShowMirResMap(self):

        resFile = self._writeMirZkAndGridResInZemax()
        writeToResMapFilePath = os.path.join(self.outputDir, "M2resMap.png")

        self.m2.showMirResMap(resFile, writeToResMapFilePath=writeToResMapFilePath)
        self.assertTrue(os.path.isfile(writeToResMapFilePath))

        os.remove(resFile)
        os.remove(writeToResMapFilePath)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
