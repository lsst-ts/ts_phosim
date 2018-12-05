import os
import numpy as np
import unittest

from lsst.ts.phosim.M1M3Sim import M1M3Sim
from lsst.ts.phosim.Utility import getModulePath


class TestM1M3Sim(unittest.TestCase):
    """ Test the M1M3Sim class."""

    def setUp(self):

        # Directory of M1M3 data
        self.mirrorDataDir = os.path.join(getModulePath(), "configData", "M1M3")

    def testFunc(self):

        # Instantiate the M1M3Sim object
        M1M3 = M1M3Sim()
        self.assertEqual(M1M3.RiInM, (2.558, 0.550))
        self.assertEqual(M1M3.RinM, (4.180, 2.508))

        M1M3.setMirrorDataDir(self.mirrorDataDir)

        forceInN = M1M3.getActForce()
        self.assertEqual(forceInN.shape, (156, 156))

        zAngleInDeg = 27.0912
        zAngleInRadian = zAngleInDeg/180*np.pi
        printthzInM = M1M3.getPrintthz(zAngleInRadian)

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM1M3Func", "M1M3printthz.txt")
        ansPrintthzInM = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(printthzInM-ansPrintthzInM)), 1e-10)

        M1M3TBulk = 0.0902
        M1M3TxGrad = -0.0894
        M1M3TyGrad = -0.1973
        M1M3TzGrad = -0.0316
        M1M3TrGrad = 0.0187
        tempCorrInUm = M1M3.getTempCorr(M1M3TBulk, M1M3TxGrad, M1M3TyGrad, M1M3TzGrad, M1M3TrGrad)

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM1M3Func", "M1M3tempCorr.txt")
        ansTempCorrInUm = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(tempCorrInUm-ansTempCorrInUm)), 6*1e-9)

        iSim = 6
        randSurfInM = M1M3.genMirSurfRandErr(zAngleInRadian, seedNum=iSim)
        
        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM1M3Func", "M1M3surfRand.txt")
        ansRandSurfInM = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(randSurfInM-ansRandSurfInM)), 1e-10)

        printthzInUm = printthzInM*1e6
        randSurfInUm = randSurfInM*1e6
        mirrorSurfInUm = printthzInUm + randSurfInUm + tempCorrInUm
        M1M3.setSurfAlongZ(mirrorSurfInUm)

        numTerms = 28
        zcInMmInZemax = M1M3.getMirrorResInMmInZemax(numTerms=numTerms)[3]

        ansFilePath = os.path.join(getModulePath(), "tests", "testData",
                                   "testM1M3Func", "sim6_M1M3zlist.txt")
        ansZcInUmInZemax = np.loadtxt(ansFilePath)
        ansZcInMmInZemax = ansZcInUmInZemax*1e-3
        self.assertLess(np.sum(np.abs(zcInMmInZemax[0:numTerms]-ansZcInMmInZemax[0:numTerms])), 1e-9)

        resFile1 = os.path.join(getModulePath(), "output", "M1res.txt")
        resFile3 = os.path.join(getModulePath(), "output", "M3res.txt")
        resFile = [resFile1, resFile3]
        M1M3.writeMirZkAndGridResInZemax(resFile=resFile, numTerms=numTerms)
        content1 = np.loadtxt(resFile1)
        content3 = np.loadtxt(resFile3)

        ansFilePath1 = os.path.join(getModulePath(), "tests", "testData",
                                    "testM1M3Func", "sim6_M1res.txt")
        ansFilePath3 = os.path.join(getModulePath(), "tests", "testData",
                                    "testM1M3Func", "sim6_M3res.txt")
        ansContent1 = np.loadtxt(ansFilePath1)
        ansContent3 = np.loadtxt(ansFilePath3)

        self.assertLess(np.sum(np.abs(content1[0,:]-ansContent1[0,:])), 1e-9)
        self.assertLess(np.sum(np.abs(content1[1:,0]-ansContent1[1:,0])), 1e-9)

        self.assertLess(np.sum(np.abs(content3[0,:]-ansContent3[0,:])), 1e-9)
        self.assertLess(np.sum(np.abs(content3[1:,0]-ansContent3[1:,0])), 1e-9)

        writeToResMapFilePath1 = os.path.join(getModulePath(), "output",
                                              "M1resMap.png")
        writeToResMapFilePath3 = os.path.join(getModulePath(), "output",
                                              "M3resMap.png")
        writeToResMapFilePath = [writeToResMapFilePath1, writeToResMapFilePath3]
        M1M3.showMirResMap(numTerms=numTerms, resFile=resFile, writeToResMapFilePath=writeToResMapFilePath)
        self.assertTrue(os.path.isfile(writeToResMapFilePath1))
        self.assertTrue(os.path.isfile(writeToResMapFilePath3))

        os.remove(resFile1)
        os.remove(resFile3)
        os.remove(writeToResMapFilePath1)
        os.remove(writeToResMapFilePath3)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
