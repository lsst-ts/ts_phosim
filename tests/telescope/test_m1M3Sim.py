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

from lsst.ts.phosim.telescope.M1M3Sim import M1M3Sim

from lsst.ts.phosim.Utility import getModulePath


class TestM1M3Sim(unittest.TestCase):
    """ Test the M1M3Sim class."""

    def setUp(self):

        self.testM3Data = os.path.join(
            getModulePath(), "tests", "testData", "testM1M3Func"
        )
        self.outputDir = os.path.join(getModulePath(), "output")

    @classmethod
    def setUpClass(cls):
        """Only do the instantiation for one time for the slow speed."""

        cls.m1m3 = M1M3Sim()

    def testInit(self):

        self.assertEqual(self.m1m3.getInnerRinM(), (2.558, 0.550))
        self.assertEqual(self.m1m3.getOuterRinM(), (4.180, 2.508))

    def testGetPrintthz(self):

        zAngleInDeg = 27.0912
        printthzInM = self._getPrintthzInM(zAngleInDeg)

        printthzMax = np.max(printthzInM)
        self.assertAlmostEqual(printthzMax * 1e7, 1.38887, places=4)
        self.assertEqual(np.where(printthzInM == printthzMax)[0][0], 3863)

        printthzMin = np.min(printthzInM)
        self.assertAlmostEqual(printthzMin * 1e7, -1.01992, places=4)
        self.assertEqual(np.where(printthzInM == printthzMin)[0][0], 3607)

        self.assertAlmostEqual(printthzInM[0] * 1e7, 1.12625, places=4)
        self.assertAlmostEqual(printthzInM[1] * 1e8, -3.40947, places=4)
        self.assertAlmostEqual(printthzInM[2] * 1e8, -2.04782, places=4)

    def _getPrintthzInM(self, zAngleInDeg):

        zAngleInRadian = np.deg2rad(zAngleInDeg)
        printthzInM = self.m1m3.getPrintthz(zAngleInRadian)

        return printthzInM

    def testGetTempCorr(self):

        tempCorrInUm = self._getTempCorrInUm()

        tempCorrMax = np.max(tempCorrInUm)
        self.assertAlmostEqual(tempCorrMax * 1e1, 9.60311, places=4)
        self.assertEqual(np.where(tempCorrInUm == tempCorrMax)[0][0], 228)

        tempCorrMin = np.min(tempCorrInUm)
        self.assertAlmostEqual(tempCorrMin * 1e1, -5.07911, places=4)
        self.assertEqual(np.where(tempCorrInUm == tempCorrMin)[0][0], 5034)

        self.assertAlmostEqual(tempCorrInUm[0] * 1e1, 4.34362, places=4)
        self.assertAlmostEqual(tempCorrInUm[1] * 1e1, 5.30053, places=4)
        self.assertAlmostEqual(tempCorrInUm[2] * 1e1, 6.49608, places=4)

    def _getTempCorrInUm(self):

        m1m3TBulk = 0.0902
        m1m3TxGrad = -0.0894
        m1m3TyGrad = -0.1973
        m1m3TzGrad = -0.0316
        m1m3TrGrad = 0.0187
        tempCorrInUm = self.m1m3.getTempCorr(
            m1m3TBulk, m1m3TxGrad, m1m3TyGrad, m1m3TzGrad, m1m3TrGrad
        )

        return tempCorrInUm

    def testGenMirSurfRandErr(self):

        iSim = 6
        zAngleInDeg = 27.0912
        randSurfInM = self._getRandSurfInM(iSim, zAngleInDeg)

        ansFilePath = os.path.join(self.testM3Data, "M1M3surfRand.txt")
        ansRandSurfInM = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(randSurfInM - ansRandSurfInM)), 1e-10)

    def _getRandSurfInM(self, iSim, zAngleInDeg):

        zAngleInRadian = np.deg2rad(zAngleInDeg)
        randSurfInM = self.m1m3.genMirSurfRandErr(zAngleInRadian, seedNum=iSim)

        return randSurfInM

    def testGetMirrorResInMmInZemax(self):

        self._setSurfAlongZ()
        zcInMmInZemax = self.m1m3.getMirrorResInMmInZemax()[3]

        ansFilePath = os.path.join(self.testM3Data, "sim6_M1M3zlist.txt")
        ansZcInUmInZemax = np.loadtxt(ansFilePath)
        ansZcInMmInZemax = ansZcInUmInZemax * 1e-3

        numTerms = self.m1m3.getNumTerms()
        delta = np.sum(np.abs(zcInMmInZemax[0:numTerms] - ansZcInMmInZemax[0:numTerms]))
        self.assertLess(delta, 1e-9)

    def _setSurfAlongZ(self):

        iSim = 6
        zAngleInDeg = 27.0912

        printthzInM = self._getPrintthzInM(zAngleInDeg)
        randSurfInM = self._getRandSurfInM(iSim, zAngleInDeg)
        tempCorrInUm = self._getTempCorrInUm()

        printthzInUm = printthzInM * 1e6
        randSurfInUm = randSurfInM * 1e6
        mirrorSurfInUm = printthzInUm + randSurfInUm + tempCorrInUm

        self.m1m3.setSurfAlongZ(mirrorSurfInUm)

    def testWriteMirZkAndGridResInZemax(self):

        resFile = self._writeMirZkAndGridResInZemax()
        resFile1, resFile3 = resFile

        content1 = np.loadtxt(resFile1)
        content3 = np.loadtxt(resFile3)

        contentShape = (41617, 4)
        self.assertEqual(content1.shape, contentShape)
        self.assertEqual(content3.shape, contentShape)

        deltaXcontent1 = content1[0, 2]
        self.assertAlmostEqual(deltaXcontent1, 42.01005, places=4)

        self.assertAlmostEqual(content1[320, 0] * 1e5, -9.97911, places=4)
        self.assertAlmostEqual(content1[320, 1] * 1e7, -1.43345, places=4)
        self.assertAlmostEqual(content1[320, 2] * 1e7, -5.70963, places=4)

        deltaXcontent3 = content3[0, 2]
        self.assertAlmostEqual(deltaXcontent3, 25.20603, places=4)

        self.assertAlmostEqual(content3[320, 0] * 1e5, -6.60547, places=4)
        self.assertAlmostEqual(content3[320, 1] * 1e8, 2.61952, places=4)
        self.assertAlmostEqual(content3[320, 2] * 1e9, -5.575, places=2)

        os.remove(resFile1)
        os.remove(resFile3)

    def _writeMirZkAndGridResInZemax(self):

        self._setSurfAlongZ()

        resFile1 = os.path.join(self.outputDir, "M1res.txt")
        resFile3 = os.path.join(self.outputDir, "M3res.txt")
        resFile = [resFile1, resFile3]
        self.m1m3.writeMirZkAndGridResInZemax(resFile=resFile)

        return resFile

    def testShowMirResMap(self):

        resFile = self._writeMirZkAndGridResInZemax()
        resFile1, resFile3 = resFile

        writeToResMapFilePath1 = os.path.join(self.outputDir, "M1resMap.png")
        writeToResMapFilePath3 = os.path.join(self.outputDir, "M3resMap.png")
        writeToResMapFilePath = [writeToResMapFilePath1, writeToResMapFilePath3]
        self.m1m3.showMirResMap(resFile, writeToResMapFilePath=writeToResMapFilePath)
        self.assertTrue(os.path.isfile(writeToResMapFilePath1))
        self.assertTrue(os.path.isfile(writeToResMapFilePath3))

        os.remove(resFile1)
        os.remove(resFile3)
        os.remove(writeToResMapFilePath1)
        os.remove(writeToResMapFilePath3)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
