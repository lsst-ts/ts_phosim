import os
import numpy as np
import unittest

from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.Utility import getModulePath


class TestOpdMetrology(unittest.TestCase):
    """ Test the OpdMetrology class."""

    def setUp(self):

        self.testDataDir = os.path.join(getModulePath(), "tests", "testData",
                                        "testOpdFunc")
        self.metr = OpdMetrology()

    def testSetWeightingRatio(self):

        wt = [1, 2]
        self.metr.setWeightingRatio(wt)
        self.assertEqual(len(self.metr.wt), len(wt))
        self.assertEqual(np.sum(self.metr.wt), 1)
        self.assertAlmostEqual(self.metr.wt[1]/self.metr.wt[0], 2)
        self.assertRaises(ValueError, self.metr.setWeightingRatio, [-1, 1])

    def testSetFieldXYinDeg(self):

        fieldXInDegree = 0.1
        fieldYInDegree = 0.2
        self.metr.setFieldXYinDeg(fieldXInDegree, fieldYInDegree)

        self.assertEqual(self.metr.fieldX, fieldXInDegree)
        self.assertEqual(self.metr.fieldY, fieldYInDegree)

    def testAddFieldXYbyDeg(self):

        fieldXInDegree = 0.1
        fieldYInDegree = 0.2
        self.metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)
        self.assertEqual(len(self.metr.fieldX), 1)

        self.metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)
        self.assertEqual(len(self.metr.fieldX), 2)

    def testSetDefaultLsstGQ(self):

        self.metr.setDefaultLsstGQ()
        self.assertEqual(len(self.metr.fieldX), 31)

    def testGetDefaultLsstWfsGQ(self):

        fieldWFSx, fieldWFSy = self.metr.getDefaultLsstWfsGQ()
        self.assertEqual(len(fieldWFSx), 4)

    def testSetDefaultComcamGQ(self):

        self.metr.setDefaultComcamGQ()
        self.assertEqual(len(self.metr.fieldX), 9)

    def testGetZkFromOpd(self):

        opdFilePath = self._getOpdFilePath()
        zk = self.metr.getZkFromOpd(opdFitsFile=opdFilePath)[0]
        
        ansOpdFileName = "sim6_iter0_opd.zer"
        ansOpdFilePath = os.path.join(self.testDataDir, ansOpdFileName)
        allOpdAns = np.loadtxt(ansOpdFilePath)
        self.assertLess(np.sum(np.abs(zk-allOpdAns[0,:])), 1e-10)

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
        xInpixel = 4000
        yInPixel = 4072
        self.metr.addFieldXYbyCamPos(sensorName, xInpixel, yInPixel,
                                     self.testDataDir)

        ansFieldXinDeg = 2000*0.2/3600
        ansFieldYinDeg = 2036*0.2/3600
        self.assertAlmostEqual((self.metr.fieldX[-1], self.metr.fieldY[-1]),
                               (ansFieldXinDeg, ansFieldYinDeg))

    def calcPSSN(self):

        pssn = self._calcPssn()
        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(pssn, allData[0,0])

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
        self.assertAlmostEqual(fwhm, allData[1,0])

    def testCalcDm5(self):

        pssn = self._calcPssn()
        dm5 = self.metr.calcDm5(pssn)

        allData = self._getMetroAllAnsData()
        self.assertAlmostEqual(dm5, allData[2,0])

    def testCalcEllip(self):

        wavelengthInUm = 0.5
        opdFilePath = self._getOpdFilePath()
        elli = self.metr.calcEllip(wavelengthInUm, opdFitsFile=opdFilePath)

        ansElliFileName = "sim6_iter0_elli.txt"
        ansElliFilePath = os.path.join(self.testDataDir, ansElliFileName)
        allElli = np.loadtxt(ansElliFilePath)
        self.assertAlmostEqual(elli, allElli[0])

    def testCalcGQvalue(self):

        self.metr.setDefaultLsstGQ()
        allData = self._getMetroAllAnsData()
        valueList = allData[0, 0:31]

        GQvalue = self.metr.calcGQvalue(valueList)
        self.assertAlmostEqual(GQvalue, allData[0,-1])


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
