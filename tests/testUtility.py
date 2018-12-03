import unittest

from lsst.ts.phosim.Utility import phosim2ZemaxCoorTrans, zemax2phosimCoorTrans


class TestUtility(unittest.TestCase):
    """Test the Utility functions."""

    def testPhosim2ZemaxCoorTrans(self):

        xInPhosim, yInPhosim, zInPhosim = 1, 2, 3
        xInZemax, yInZemax, zInZemax = phosim2ZemaxCoorTrans(
                                        xInPhosim, yInPhosim, zInPhosim)
        self.assertEqual((xInZemax, yInZemax, zInZemax),
                         (-xInPhosim, yInPhosim, -zInPhosim))

    def testZemax2phosimCoorTrans(self):

        xInZemax, yInZemax, zInZemax = 1, 2, 3
        xInPhosim, yInPhosim, zInPhosim = zemax2phosimCoorTrans(
                                            xInZemax, yInZemax, zInZemax)
        self.assertEqual((xInPhosim, yInPhosim, zInPhosim),
                         (-xInZemax, yInZemax, -zInZemax))


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
