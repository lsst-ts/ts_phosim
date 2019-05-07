import os
import unittest
import numpy as np

from lsst.ts.phosim.MirrorSim import MirrorSim
from lsst.ts.phosim.Utility import getConfigDir


class TestMirrorSim(unittest.TestCase):
    """ Test the MirrorSim class."""

    def setUp(self):

        self.innerRinM = 0.9
        self.outerRinM = 1.710

        configDir = os.path.join(getConfigDir(), "M2")
        self.mirror = MirrorSim(self.innerRinM, self.outerRinM, configDir)
        self.mirror.config(numTerms=28,
                           actForceFileName="M2_1um_force.yaml",
                           lutFileName="")

    def testInit(self):

        self.assertEqual(self.mirror.getInnerRinM(), self.innerRinM)
        self.assertEqual(self.mirror.getOuterRinM(), self.outerRinM)

    def testGetNumTerms(self):

        self.assertEqual(self.mirror.getNumTerms(), 28)

    def testSetAndGetSurfAlongZ(self):

        surfAlongZinUm = np.random.rand(3, 3)
        self.mirror.setSurfAlongZ(surfAlongZinUm)

        delta = np.sum(np.abs(self.mirror.getSurfAlongZ() - surfAlongZinUm))
        self.assertEqual(delta, 0)

    def testGetLUTforce(self):

        m1m3DataDir = os.path.join(getConfigDir(), "M1M3")
        mirror = MirrorSim(self.innerRinM, self.outerRinM, m1m3DataDir)
        mirror.config(numTerms=28, lutFileName="M1M3_LUT.yaml")

        zangleInDeg = 1.5
        lutForce = mirror.getLUTforce(zangleInDeg)

        # Skip the row 1 as the ruler
        oriLutForce = mirror._lutFile.getMatContent()[1:, :]
        ansLutForce = (oriLutForce[:, 1]+oriLutForce[:, 2])/2

        self.assertLess(np.sum(np.abs(lutForce-ansLutForce)), 1e-10)

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

    def testShowMirResMap(self):

        self.assertRaises(NotImplementedError,
                          self.mirror.showMirResMap, "")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
