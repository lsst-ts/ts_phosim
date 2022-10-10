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
import unittest

import numpy as np
from lsst.ts.phosim.utils.ConvertZernikesToPsfWidth import (
    convertZernikesToPsfWidth,
    getPsfGradPerZernike,
)
from lsst.ts.phosim.utils.Utility import getModulePath
from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.Utility import CamType


class TestConvertZernikesToPsfWidth(unittest.TestCase):
    """Test the convertZernikesToPsfWidth function."""

    def setUp(self):
        """Load saved conversion factors."""
        self.testDataDir = os.path.join(
            getModulePath(), "tests", "testData", "psfGradientsPerZernike"
        )

    def testLsstCam(self):
        """Test that the LsstCam values match the expected values."""
        # LsstCam should be selected by default
        conversion_factors = getPsfGradPerZernike(37)
        expected_factors = np.genfromtxt(os.path.join(self.testDataDir, "lsstcam.txt"))
        self.assertTrue(np.allclose(conversion_factors, expected_factors, atol=1e-3))

    def testAuxTel(self):
        """Test that the AuxTel values match the expected values."""
        # Setup the AuxTel instrument
        # Note the donut dimension shouldn't matter for this computation
        inst = Instrument()
        inst.configFromFile(160, CamType.AuxTel)

        # Calculate and compare conversion factors
        conversion_factors = getPsfGradPerZernike(37, inst)
        expected_factors = np.genfromtxt(os.path.join(self.testDataDir, "auxtel.txt"))
        self.assertTrue(np.allclose(conversion_factors, expected_factors, atol=1e-3))

    def testAllCamTypes(self):
        """Test that we can convert Zernikes for all cameras."""
        dummy_zernikes = np.ones(24)
        for cam in CamType:
            with self.subTest(cam=cam):
                # Create the instrument
                inst = Instrument()
                inst.configFromFile(160, cam)

                # Convert zernikes
                converted_zernikes = convertZernikesToPsfWidth(dummy_zernikes, inst)

                self.assertTrue(np.all(converted_zernikes > 0))


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
