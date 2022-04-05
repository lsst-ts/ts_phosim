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

from lsst.ts.phosim.telescope.CamSim import CamSim

from lsst.ts.phosim.utils.Utility import getModulePath, CamDistType


class TestCamSim(unittest.TestCase):
    """Test the CamSim class."""

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
            zAngleInRad, CamDistType.L1S1zer
        )

        dataFilePath = os.path.join(
            getModulePath(), "tests", "testData", "testOpdFunc", "sim6_iter0_pert.cmd"
        )
        distData = np.loadtxt(dataFilePath, skiprows=88, usecols=(1, 2, 3))
        idx = distData[:, 0] == 3
        absDiff = np.sum(np.abs(distortionInMn - distData[idx, -1]))

        self.assertTrue(absDiff < 1e-10)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
