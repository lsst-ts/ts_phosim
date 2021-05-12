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

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getModulePath


class TestSkySim(unittest.TestCase):
    """ Test the SkySim class."""

    def setUp(self):

        self.skySim = SkySim()

    def testGetStarId(self):

        starId = self.skySim.getStarId()
        self.assertEqual(len(starId), 0)
        self.assertTrue(isinstance(starId, np.ndarray))

    def testGetRaDecInDeg(self):

        ra, decl = self.skySim.getRaDecInDeg()

        self.assertEqual(len(ra), 0)
        self.assertTrue(isinstance(ra, np.ndarray))

        self.assertEqual(len(decl), 0)
        self.assertTrue(isinstance(decl, np.ndarray))

    def testGetStarMag(self):

        mag = self.skySim.getStarMag()
        self.assertEqual(len(mag), 0)
        self.assertTrue(isinstance(mag, np.ndarray))

    def testAddStarByRaDecInDeg(self):

        self.skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.assertEqual(len(self.skySim.getStarId()), 1)

        self.skySim.addStarByRaDecInDeg(2, 2.1, 3, 4)
        starId = self.skySim.getStarId()
        self.assertEqual(len(starId), 2)
        self.assertEqual(starId[0], 1)
        self.assertEqual(starId[1], 2)

        # Try to add the same star Id again
        self.skySim.addStarByRaDecInDeg(2, 2.1, 3, 4)
        self.assertEqual(len(self.skySim.getStarId()), 2)

    def testResetSky(self):

        self.skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.skySim.resetSky()
        self.assertEqual(len(self.skySim.getStarId()), 0)

    def testSetStarRaDecInDeg(self):

        self.skySim.setStarRaDecInDeg(
            np.array([0]), np.array([1]), np.array([2]), np.array([3])
        )
        self.assertEqual(len(self.skySim.getStarId()), 1)

    def testAddStarByFile(self):

        self._addStarByFile("wfsStar.txt")

        self.assertEqual(len(self.skySim.getStarId()), 8)

        ra, decl = self.skySim.getRaDecInDeg()
        self.assertEqual(ra[2], -1.176)
        self.assertEqual(decl[2], 1.196)

        self.assertEqual(self.skySim.getStarMag()[2], 17.0)

    def _addStarByFile(self, skyFileName):

        skyFile = os.path.join(getModulePath(), "tests", "testData", "sky", skyFileName)
        self.skySim.addStarByFile(skyFile)

    def testAddStarByFileWithSglStar(self):

        self._addStarByFile("wfsSglStar.txt")

        self.assertEqual(len(self.skySim.getStarId()), 1)

        ra, decl = self.skySim.getRaDecInDeg()
        self.assertEqual(ra[0], 1.196)
        self.assertEqual(decl[0], 1.176)

        self.assertEqual(self.skySim.getStarMag()[0], 17.0)

    def testExportSkyToFile(self):

        self._addStarByFile("wfsStar.txt")
        outputFilePath = os.path.join(getModulePath(), "output", "testSkyOutput.txt")

        self.skySim.exportSkyToFile(outputFilePath)
        self.assertTrue(os.path.isfile(outputFilePath))
        os.remove(outputFilePath)

    def testAddStarByChipPos(self):

        self._setObservationMetaData()

        # Add the star
        sensorName = "R22_S11"
        starId = 0
        xInpixelInCam = 2000
        yInPixelInCam = 2036
        starMag = 17
        self.skySim.addStarByChipPos(
            sensorName, starId, xInpixelInCam, yInPixelInCam, starMag
        )

        # Test the result
        ra, decl = self.skySim.getRaDecInDeg()
        self.assertAlmostEqual(ra[0], 0.00063889)
        self.assertAlmostEqual(decl[0], -8.333e-05)

    def _setObservationMetaData(self):

        ra = 0
        decl = 0
        rotSkyPos = 0
        self.skySim.setObservationMetaData(ra, decl, rotSkyPos)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
