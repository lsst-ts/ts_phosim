import os
import unittest
import numpy as np

from lsst.ts.phosim.M2Sim import M2Sim
from lsst.ts.phosim.Utility import getModulePath


class TestM2Sim(unittest.TestCase):
    """ Test the M2Sim class."""

    def setUp(self):

        mirrorDataDir = os.path.join(getModulePath(), "configData", "M2")
        self.M2 = M2Sim(mirrorDataDir=mirrorDataDir)

    def testInit(self):

        self.assertEqual(self.M2.getInnerRinM(), 0.9)
        self.assertEqual(self.M2.getOuterRinM(), 1.71)

    def testGetActForce(self):

        forceInN = self.M2.getActForce()
        self.assertEqual(forceInN.shape, (156, 156))

    def testGetPrintthz(self):

        zAngleInDeg = 27.0912
        zAngleInRadian = np.deg2rad(zAngleInDeg)
        printthzInUm = self.M2.getPrintthz(zAngleInRadian)

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM2Func", "M2printthz.txt")
        ansPrintthzInUm = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(printthzInUm-ansPrintthzInUm)), 1e-10)

    def testGetTempCorr(self):

        M2TzGrad = -0.0675
        M2TrGrad = -0.1416
        tempCorrInUm = self.M2.getTempCorr(M2TzGrad, M2TrGrad)

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM2Func", "M2tempCorr.txt")
        ansTempCorrInUm = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(tempCorrInUm-ansTempCorrInUm)), 1e-10)

    def testGetMirrorResInMmInZemax(self):

        numTerms = 28
        self._setSurfAlongZ()
        zcInMmInZemax = self.M2.getMirrorResInMmInZemax(numTerms=numTerms)[3]

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM2Func", "sim6_M2zlist.txt")
        ansZcInUmInZemax = np.loadtxt(ansFilePath)
        ansZcInMmInZemax = ansZcInUmInZemax*1e-3

        delta = np.sum(np.abs(zcInMmInZemax[0:numTerms] -
                              ansZcInMmInZemax[0:numTerms]))
        self.assertLess(delta, 1e-9)

    def _setSurfAlongZ(self):

        zAngleInDeg = 27.0912
        zAngleInRadian = np.deg2rad(zAngleInDeg)
        printthzInUm = self.M2.getPrintthz(zAngleInRadian)

        M2TzGrad = -0.0675
        M2TrGrad = -0.1416
        tempCorrInUm = self.M2.getTempCorr(M2TzGrad, M2TrGrad)

        mirrorSurfInUm = printthzInUm + tempCorrInUm
        self.M2.setSurfAlongZ(mirrorSurfInUm)

    def testWriteMirZkAndGridResInZemax(self):

        numTerms = 28
        self._setSurfAlongZ()
        resFile = os.path.join(getModulePath(), "output", "M2res.txt")
        self.M2.writeMirZkAndGridResInZemax(resFile=resFile, numTerms=numTerms)
        content = np.loadtxt(resFile)

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM2Func", "sim6_M2res.txt")
        ansContent = np.loadtxt(ansFilePath)

        self.assertLess(np.sum(np.abs(content[0, :]-ansContent[0, :])), 1e-9)
        self.assertLess(np.sum(np.abs(content[1:, 0]-ansContent[1:, 0])), 1e-9)

        os.remove(resFile)

    def testShowMirResMap(self):

        numTerms = 28
        self._setSurfAlongZ()
        resFile = os.path.join(getModulePath(), "output", "M2res.txt")
        self.M2.writeMirZkAndGridResInZemax(resFile=resFile, numTerms=numTerms)

        writeToResMapFilePath = os.path.join(getModulePath(), "output",
                                             "M2resMap.png")
        self.M2.showMirResMap(numTerms=numTerms, resFile=resFile,
                              writeToResMapFilePath=writeToResMapFilePath)
        self.assertTrue(os.path.isfile(writeToResMapFilePath))

        os.remove(resFile)
        os.remove(writeToResMapFilePath)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
