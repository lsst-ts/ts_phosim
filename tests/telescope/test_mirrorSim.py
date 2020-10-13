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

from lsst.ts.phosim.telescope.MirrorSim import MirrorSim

from lsst.ts.phosim.Utility import getConfigDir


class TestMirrorSim(unittest.TestCase):
    """ Test the MirrorSim class."""

    @classmethod
    def setUpClass(cls):
        """Only do the instantiation for one time for the slow speed."""

        cls.innerRinM = 0.9
        cls.outerRinM = 1.710

        configDir = os.path.join(getConfigDir(), "M2")
        cls.mirror = MirrorSim(cls.innerRinM, cls.outerRinM, configDir)
        cls.mirror.config(numTerms=28, lutFileName="")

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
        ansLutForce = (oriLutForce[:, 1] + oriLutForce[:, 2]) / 2

        self.assertLess(np.sum(np.abs(lutForce - ansLutForce)), 1e-10)

    def testGetPrintthz(self):

        self.assertRaises(NotImplementedError, self.mirror.getPrintthz, 0)

    def testGetTempCorr(self):

        self.assertRaises(NotImplementedError, self.mirror.getTempCorr)

    def testGetMirrorResInMmInZemax(self):

        self.assertRaises(NotImplementedError, self.mirror.getMirrorResInMmInZemax)

    def testWriteMirZkAndGridResInZemax(self):

        self.assertRaises(NotImplementedError, self.mirror.writeMirZkAndGridResInZemax)

    def testShowMirResMap(self):

        self.assertRaises(NotImplementedError, self.mirror.showMirResMap, "")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
