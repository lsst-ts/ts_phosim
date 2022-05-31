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

import unittest

import numpy as np
from lsst.ts.phosim.utils.ConvertZernikesToPsfWidth import (
    _load_conversion_factors,
    convertZernikesToPsfWidth,
)
from lsst.ts.wep.Utility import CamType


class TestConvertZernikesToPsfWidth(unittest.TestCase):
    """Test the convertZernikesToPsfWidth function."""

    def setUp(self):
        """Load conversion factors, and sets of valid and invalid cameras."""
        self.conversion_factors = _load_conversion_factors()
        self.valid_cameras = set(self.conversion_factors)
        self.invalid_cameras = set(CamType) - self.valid_cameras

    def testOnes(self):
        """Test that converting array of 1s returns the conversion factors."""
        for camType in self.valid_cameras:
            with self.subTest(camType=camType):
                conversion_factors = self.conversion_factors[camType]
                zks = np.ones(len(conversion_factors))
                zks_arcsecs = convertZernikesToPsfWidth(zks, camType)
                self.assertTrue(np.allclose(zks_arcsecs, conversion_factors))

    def testTooManyZernikes(self):
        """Test that ValueError is raised when too many zernikes are passed."""
        for camType in self.valid_cameras:
            with self.subTest(camType=camType):
                conversion_factors = self.conversion_factors[camType]
                zks = np.ones(len(conversion_factors) + 1)
                with self.assertRaises(ValueError):
                    convertZernikesToPsfWidth(zks, camType)

    def testInvalidCameras(self):
        """Test that supplying unsupported camType raises a KeyError."""
        for camType in self.invalid_cameras:
            with self.subTest(camType=camType):
                with self.assertRaises(KeyError):
                    convertZernikesToPsfWidth(np.ones(1), camType)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
