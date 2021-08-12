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

from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.Utility import getModulePath


class TestOpdMetrology(unittest.TestCase):
    """ Test the OpdMetrology class."""

    def setUp(self):

        self.testDataDir = os.path.join(
            getModulePath(), "tests", "testData", "testOpdFunc"
        )
        self.metr = OpdMetrology()

    def testGetFieldXY(self):

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(len(fieldX), 0)
        self.assertEqual(len(fieldY), 0)
        self.assertTrue(isinstance(fieldX, np.ndarray))
        self.assertTrue(isinstance(fieldY, np.ndarray))

    def testGetWeightingRatio(self):

        wt = self.metr.getWeightingRatio()
        self.assertEqual(len(wt), 0)
        self.assertTrue(isinstance(wt, np.ndarray))

    def testSetWeightingRatio(self):

        wt = [1, 2]
        self.metr.setWeightingRatio(wt)

        wtInMetr = self.metr.getWeightingRatio()
        self.assertEqual(len(wtInMetr), len(wt))
        self.assertEqual(np.sum(wtInMetr), 1)
        self.assertAlmostEqual(wtInMetr[1] / wtInMetr[0], 2)
        self.assertRaises(ValueError, self.metr.setWeightingRatio, [-1, 1])

    def testSetFieldXYinDeg(self):

        fieldXInDegree = 0.1
        fieldYInDegree = 0.2
        self.metr.setFieldXYinDeg(fieldXInDegree, fieldYInDegree)

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(fieldX, fieldXInDegree)
        self.assertEqual(fieldY, fieldYInDegree)

    def testAddFieldXYbyDeg(self):

        fieldXInDegree = 0.1
        fieldYInDegree = 0.2
        self.metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(len(fieldX), 1)

        self.metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(len(fieldX), 2)

    def testSetWgtAndFieldXyOfGQLsst(self):

        self.metr.setWgtAndFieldXyOfGQ("lsst")

        fieldXAns, fieldYAns, wgtAns = self._calcFieldXyAndWgtLsst()

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(len(fieldX), 31)
        self.assertLess(np.sum(np.abs(fieldX - fieldXAns)), 1e-10)
        self.assertLess(np.sum(np.abs(fieldY - fieldYAns)), 1e-10)

        wgt = self.metr.getWeightingRatio()
        self.assertEqual(len(wgt), 31)
        self.assertLess(np.sum(np.abs(wgt - wgtAns)), 1e-10)

    def _calcFieldXyAndWgtLsst(self):

        # The distance of point xi (used in Gaussian quadrature plane) to the
        # origin
        # This value is in [-1.75, 1.75]
        armLen = [0.379, 0.841, 1.237, 1.535, 1.708]

        # Weighting of point xi (used in Gaussian quadrature plane) for each
        # ring
        armW = [0.2369, 0.4786, 0.5689, 0.4786, 0.2369]

        # Number of points on each ring
        nArm = 6

        # Get the weighting for all field points (31 for lsst camera)
        # Consider the first element is center (0)
        wgt = np.concatenate([np.zeros(1), np.kron(armW, np.ones(nArm))])

        # Generate the fields point x, y coordinates
        pointAngle = np.arange(nArm) * (2 * np.pi) / nArm
        fieldX = np.concatenate([np.zeros(1), np.kron(armLen, np.cos(pointAngle))])
        fieldY = np.concatenate([np.zeros(1), np.kron(armLen, np.sin(pointAngle))])

        return fieldX, fieldY, wgt / np.sum(wgt)

    def testSetWgtAndFieldXyOfGQComCam(self):

        self.metr.setWgtAndFieldXyOfGQ("comcam")

        fieldXAns, fieldYAns, wgtAns = self._calcFieldXyAndWgtComCam()

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(len(fieldX), 9)
        self.assertLess(np.sum(np.abs(fieldX - fieldXAns)), 1e-10)
        self.assertLess(np.sum(np.abs(fieldY - fieldYAns)), 1e-10)

        wgt = self.metr.getWeightingRatio()
        self.assertEqual(len(wgt), 9)
        self.assertLess(np.sum(np.abs(wgt - wgtAns)), 1e-10)

    def _calcFieldXyAndWgtComCam(self):

        # ComCam is the cetral raft of LSST cam, which is composed of 3 x 3
        # CCDs.
        nRow = 3
        nCol = 3

        # Number of field points
        nField = nRow * nCol

        # Get the weighting for all field points (9 for comcam)
        wgt = np.ones(nField)

        # Distance to raft center in degree along x/y direction and the
        # related relative position
        sensorD = 0.2347
        coorComcam = sensorD * np.array([-1, 0, 1])

        # Generate the fields point x, y coordinates
        fieldX = np.kron(coorComcam, np.ones(nRow))
        fieldY = np.kron(np.ones(nCol), coorComcam)

        return fieldX, fieldY, wgt / np.sum(wgt)

    def testSetWgtAndFieldXyOfGQLsstFam(self):

        self.metr.setWgtAndFieldXyOfGQ("lsstfam")

        fieldX, fieldY = self.metr.getFieldXY()
        self.assertEqual(len(fieldX), 189)

        wgt = self.metr.getWeightingRatio()
        self.assertEqual(len(wgt), 189)

    def testSetWgtAndFieldXyOfGQErr(self):

        self.assertRaises(
            RuntimeError, self.metr.setWgtAndFieldXyOfGQ, "NoThisInstName"
        )

    def testGetDefaultLsstWfsGQ(self):

        fieldWFSx, fieldWFSy = self.metr.getDefaultLsstWfsGQ()
        self.assertEqual(len(fieldWFSx), 4)

    def testGetZkFromOpd(self):

        opdFilePath = self._getOpdFilePath()
        zk = self.metr.getZkFromOpd(opdFitsFile=opdFilePath)[0]

        ansOpdFileName = "sim6_iter0_opd.zer"
        ansOpdFilePath = os.path.join(self.testDataDir, ansOpdFileName)
        allOpdAns = np.loadtxt(ansOpdFilePath)
        self.assertLess(np.sum(np.abs(zk - allOpdAns[0, :])), 1e-10)

    def _getOpdFilePath(self):

        opdFileName = "sim6_iter0_opd0.fits.gz"
        opdFilePath = os.path.join(self.testDataDir, opdFileName)

        return opdFilePath

    def testRmPTTfromOPD(self):

        opdFilePath = self._getOpdFilePath()
        opdRmPTT, opdx, opdy = self.metr.rmPTTfromOPD(opdFitsFile=opdFilePath)

        zkRmPTT = self.metr.getZkFromOpd(opdMap=opdRmPTT)[0]
        self.assertLess(np.sum(np.abs(zkRmPTT[0:3])), 5e-2)

    def testAddFieldXYbyCamPos(self):

        sensorName = "R22_S11"
        xInpixel = 4004
        yInPixel = 4096
        self.metr.addFieldXYbyCamPos(sensorName, xInpixel, yInPixel, self.testDataDir)

        ansFieldXinDeg = 2002 * 0.2 / 3600
        ansFieldYinDeg = 2048 * 0.2 / 3600
        fieldX, fieldY = self.metr.getFieldXY()
        self.assertAlmostEqual(
            (fieldX[-1], fieldY[-1]), (ansFieldXinDeg, ansFieldYinDeg)
        )

    def calcPSSN(self):

        pssn = self._calcPssn()
        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(pssn, allData[0, 0])

    def _calcPssn(self):

        wavelengthInUm = 0.5
        opdFilePath = self._getOpdFilePath()
        pssn = self.metr.calcPSSN(wavelengthInUm, opdFitsFile=opdFilePath)

        return pssn

    def _getMetroAllAnsData(self):

        ansAllDataFileName = "sim6_iter0_PSSN.txt"
        ansAllDataFilePath = os.path.join(self.testDataDir, ansAllDataFileName)
        allData = np.loadtxt(ansAllDataFilePath)

        return allData

    def testCalcFWHMeff(self):

        pssn = self._calcPssn()
        fwhm = self.metr.calcFWHMeff(pssn)

        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(fwhm, allData[1, 0])

    def testCalcDm5(self):

        pssn = self._calcPssn()
        dm5 = self.metr.calcDm5(pssn)

        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(dm5, allData[2, 0])

    def testCalcEllip(self):

        wavelengthInUm = 0.5
        opdFilePath = self._getOpdFilePath()
        elli = self.metr.calcEllip(wavelengthInUm, opdFitsFile=opdFilePath)

        ansElliFileName = "sim6_iter0_elli.txt"
        ansElliFilePath = os.path.join(self.testDataDir, ansElliFileName)
        allElli = np.loadtxt(ansElliFilePath)
        self.assertAlmostEqual(elli, allElli[0])

    def testCalcGQvalue(self):

        self.metr.setWgtAndFieldXyOfGQ("lsst")
        allData = self._getMetroAllAnsData()
        valueList = allData[0, 0:31]

        GQvalue = self.metr.calcGQvalue(valueList)
        self.assertAlmostEqual(GQvalue, allData[0, -1])


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
