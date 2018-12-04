import os
import unittest
import numpy as np

from lsst.ts.phosim.MirrorSim import MirrorSim
from lsst.ts.phosim.Utility import getModulePath


class TestMirrorSim(unittest.TestCase):
    """ Test the MirrorSim class."""

    def setUp(self):

        # Directory of M2 data
        self.mirrorDataDir = os.path.join(getModulePath(), "configData", "M2")

        # Directory of M1M3 data
        self.M1M3DataDir = os.path.join(getModulePath(), "configData", "M1M3")

    def testFunc(self):

        # Inner raidus
        innerRinM = 0.9

        # Outer raidus
        outerRinM = 1.710

        # Instantiate the MirrorSim object
        mirror = MirrorSim(innerRinM, outerRinM)
        self.assertEqual(mirror.RiInM, innerRinM)
        self.assertEqual(mirror.RinM, outerRinM)

        mirror.setMirrorDataDir(self.mirrorDataDir)
        self.assertEqual(mirror.mirrorDataDir, self.mirrorDataDir)

        dataFileName = "M2_GT_FEA.txt"
        data = mirror.getMirrorData(dataFileName, skiprows=1)
        self.assertEqual(data.shape, (9084, 6))

        surfAlongZ = np.random.rand(3, 4)
        mirror.setSurfAlongZ(surfAlongZ)
        self.assertEqual(np.sum(np.abs(mirror.surf-surfAlongZ)), 0)

        zangleInDeg = 1.5
        LUTfileName = "M1M3_LUT.txt"
        mirror.setMirrorDataDir(self.M1M3DataDir)
        LUTforce = mirror.getLUTforce(zangleInDeg, LUTfileName)

        oriLutForce = mirror.getMirrorData(LUTfileName, skiprows=1)
        ansLutForce = (oriLutForce[:,1]+oriLutForce[:,2])/2
        self.assertLess(np.sum(np.abs(LUTforce-ansLutForce)), 1e-10)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
