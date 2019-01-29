import unittest

from lsst.ts.wep.Utility import FilterType
from lsst.ts.phosim.Utility import opt2ZemaxCoorTrans, zemax2optCoorTrans, \
                                   mapSurfNameToEnum, SurfaceType, getWfChips, \
                                   createObservation, getOpsimObservation


class TestUtility(unittest.TestCase):
    """Test the Utility functions."""

    def testGetWfChips(self):

        wfChips = getWfChips()
        self.assertEqual(len(wfChips), 4)
        self.assertEqual(len(wfChips[0]), 2)

    def testCreateObservationFullParameters(self):

        obsId = 100
        filterType = FilterType.U
        boresight = (10, 20)
        zAngleInDeg = 10.0
        rotAngInDeg = 11.0
        mjd = 4000.0

        obs = createObservation(obsId=obsId, filterType=filterType,
                            boresight=boresight, zAngleInDeg=zAngleInDeg,
                            rotAngInDeg=rotAngInDeg, mjd=mjd)

        self.assertEqual(obs.OpsimMetaData["obsHistID"], obsId)
        self.assertEqual(obs.bandpass, filterType.toString())
        self.assertEqual(obs.pointingRA, boresight[0])
        self.assertEqual(obs.pointingDec, boresight[1])
        self.assertEqual(90 - obs.OpsimMetaData["altitude"], zAngleInDeg)
        self.assertEqual(obs.rotSkyPos, rotAngInDeg)
        self.assertAlmostEqual(obs.mjd.TAI, mjd)

    def testCreateObservationMissingParameters(self):

        obsId = 100
        obs = createObservation(obsId=obsId)

        self.assertEqual(obs.OpsimMetaData["obsHistID"], obsId)
        self.assertEqual(obs.rotSkyPos, 0)

    def testCreateObservationBadZenith(self):

        badZenith = -45
        with self.assertRaises(ValueError):
            createObservation(zAngleInDeg=badZenith)

    def testGetOpsimObservation(self):

        obs = getOpsimObservation(0, 0)
        self.assertEqual(obs.bandpass, 'y')

    def testGetOpsimObservationBadTract(self):

        badTract = -1
        goodTarget = 10
        with self.assertRaises(ValueError):
            getOpsimObservation(badTract, goodTarget)

    def testGetOpsimObservationBadTarget(self):

        goodTract = 0
        badTarget = 50
        with self.assertRaises(ValueError):
            getOpsimObservation(goodTract, badTarget)

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
