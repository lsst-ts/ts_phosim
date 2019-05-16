import os
import numpy as np
import unittest

from lsst.ts.phosim.telescope.CamSim import CamSim

from lsst.ts.phosim.Utility import getModulePath, CamDistType


class TestCamSim(unittest.TestCase):
    """ Test the CamSim class."""

    def setUp(self):

        self.camSim = CamSim()

    def testSetRotAngInRad(self):

        rotAngInRad = 1.0
        self.camSim.setRotAngInRad(rotAngInRad)
        self.assertEqual(self.camSim.camRotInRad, rotAngInRad)

        self.assertRaises(ValueError, self.camSim.setRotAngInRad, -4)
        self.assertRaises(ValueError, self.camSim.setRotAngInRad, 4)

    def testSetRotAngInDeg(self):

        rotAngInDeg = 30.0
        self.camSim.setRotAngInDeg(rotAngInDeg)
        self.assertEqual(self.camSim.camRotInRad, np.deg2rad(rotAngInDeg))

        self.assertRaises(ValueError, self.camSim.setRotAngInDeg, -91)
        self.assertRaises(ValueError, self.camSim.setRotAngInDeg, 91)

    def testSetBodyTempInDegC(self):

        tempInDegC = 10.0
        self.camSim.setBodyTempInDegC(tempInDegC)
        self.assertEqual(self.camSim.camTBinDegC, tempInDegC)

        self.assertRaises(ValueError, self.camSim.setBodyTempInDegC, 1.9)
        self.assertRaises(ValueError, self.camSim.setBodyTempInDegC, 16.1)

    def testGetCamDistortionInMm(self):

        self.camSim.setBodyTempInDegC(6.5650)
        self.camSim.setRotAngInRad(-1.2323)

        zenithAngleInDeg = 27.0912
        zAngleInRad = np.deg2rad(zenithAngleInDeg)
        distortionInMn = self.camSim.getCamDistortionInMm(
            zAngleInRad, CamDistType.L1S1zer)

        dataFilePath = os.path.join(getModulePath(), "tests", "testData",
                                    "testOpdFunc", "sim6_iter0_pert.cmd")
        distData = np.loadtxt(dataFilePath, skiprows=88, usecols=(1, 2, 3))
        idx = (distData[:, 0] == 3)
        absDiff = np.sum(np.abs(distortionInMn - distData[idx, -1]))

        self.assertTrue(absDiff < 1e-10)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
