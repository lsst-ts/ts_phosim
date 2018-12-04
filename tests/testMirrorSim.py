import os
import unittest
import numpy as np

from lsst.ts.phosim.MirrorSim import MirrorSim
from lsst.ts.phosim.Utility import getModulePath


class TestMirrorSim(unittest.TestCase):
    """ Test the MirrorSim class."""

    def setUp(self):

        self.innerRinM = 0.9
        self.outerRinM = 1.710
        self.mirror = MirrorSim(self.innerRinM, self.outerRinM)

    def testInit(self):

        self.assertEqual(self.mirror.getInnerRinM(), self.innerRinM)
        self.assertEqual(self.mirror.getOuterRinM(), self.outerRinM)

    def testSetMirrorDataDir(self):

        mirrorDataDir = "MirrorDataDir"
        self.mirror.setMirrorDataDir(mirrorDataDir)
        self.assertEqual(self.mirror.mirrorDataDir, mirrorDataDir)

    def testGetMirrorData(self):

        M2DataDir = os.path.join(getModulePath(), "configData", "M2")
        self.mirror.setMirrorDataDir(M2DataDir)

        dataFileName = "M2_GT_FEA.txt"
        data = self.mirror.getMirrorData(dataFileName, skiprows=1)
        self.assertEqual(data.shape, (9084, 6))

    def testSetAndGetSurfAlongZ(self):

        surfAlongZinUm = np.random.rand(3, 3)
        self.mirror.setSurfAlongZ(surfAlongZinUm)

        delta = np.sum(np.abs(self.mirror.getSurfAlongZ() - surfAlongZinUm))
        self.assertEqual(delta, 0)

    def testGetLUTforce(self):

        M1M3DataDir = os.path.join(getModulePath(), "configData", "M1M3")

        zangleInDeg = 1.5
        LUTfileName = "M1M3_LUT.txt"
        self.mirror.setMirrorDataDir(M1M3DataDir)
        LUTforce = self.mirror.getLUTforce(zangleInDeg, LUTfileName)

        oriLutForce = self.mirror.getMirrorData(LUTfileName, skiprows=1)
        ansLutForce = (oriLutForce[:,1]+oriLutForce[:,2])/2
        self.assertLess(np.sum(np.abs(LUTforce-ansLutForce)), 1e-10)

    def testGetPrintthz(self):

        self.assertRaises(NotImplementedError, self.mirror.getPrintthz, 0)

    def testGetTempCorr(self):

        self.assertRaises(NotImplementedError, self.mirror.getTempCorr)

    def testGetMirrorResInMmInZemax(self):

        self.assertRaises(NotImplementedError,
                          self.mirror.getMirrorResInMmInZemax)

    def testWriteMirZkAndGridResInZemax(self):

        self.assertRaises(NotImplementedError,
                          self.mirror.writeMirZkAndGridResInZemax)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
