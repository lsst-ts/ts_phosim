import unittest

from lsst.ts.phosim.Utility import opt2ZemaxCoorTrans, zemax2optCoorTrans, \
    mapSurfNameToEnum, SurfaceType


class TestUtility(unittest.TestCase):
    """Test the Utility functions."""

    def testPhosim2ZemaxCoorTrans(self):

        xInPhosim, yInPhosim, zInPhosim = 1, 2, 3
        xInZemax, yInZemax, zInZemax = opt2ZemaxCoorTrans(
            xInPhosim, yInPhosim, zInPhosim)
        self.assertEqual((xInZemax, yInZemax, zInZemax),
                         (-xInPhosim, yInPhosim, -zInPhosim))

    def testZemax2phosimCoorTrans(self):

        xInZemax, yInZemax, zInZemax = 1, 2, 3
        xInPhosim, yInPhosim, zInPhosim = zemax2optCoorTrans(
            xInZemax, yInZemax, zInZemax)
        self.assertEqual((xInPhosim, yInPhosim, zInPhosim),
                         (-xInZemax, yInZemax, -zInZemax))

    def testMapSurfNameToEnum(self):

        surfName = "L2F"
        surfaceType = mapSurfNameToEnum(surfName)
        self.assertEqual(surfaceType, SurfaceType.L2F)
        self.assertRaises(ValueError, mapSurfNameToEnum, "L123F")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
