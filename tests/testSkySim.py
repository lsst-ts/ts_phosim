import os
import numpy as np
import unittest

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getModulePath, FilterType


class TestSkySim(unittest.TestCase):
    """ Test the SkySim class."""

    def setUp(self):

        # Directory to the focal plane file
        folderPath2FocalPlane = os.path.join(getModulePath(), "tests",
                                             "testData", "testOpdFunc")
        self.skySim = SkySim()
        self.skySim.setFolderPath2FocalPlane(folderPath2FocalPlane)

    def testConfigDbInfo(self):

        self._setDbInfo()
        self.assertEqual(self.skySim.dbInfo["host"], "localhost:51433")

    def _setDbInfo(self):

        databaseHost = "localhost:51433"
        databaseUser = "LSST-2"
        databasePassword = "L$$TUser"
        databaseName = "LSSTCATSIM"

        self.skySim.configDbInfo(databaseHost, databaseUser,
                                 databasePassword, databaseName)

    def testAddStarByRaDecInDeg(self):

        self.skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.assertEqual(len(self.skySim.starId), 1)

        self.skySim.addStarByRaDecInDeg(2, 2.1, 3, 4)
        self.assertEqual(len(self.skySim.starId), 2)
        self.assertEqual(self.skySim.starId[0], 1)
        self.assertEqual(self.skySim.starId[1], 2)

        # Try to add the same star Id again
        self.skySim.addStarByRaDecInDeg(2, 2.1, 3, 4)
        self.assertEqual(len(self.skySim.starId), 2)

    def testResetSky(self):

        self.skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.skySim.resetSky()
        self.assertEqual(len(self.skySim.starId), 0)

    def testSetStarRaDecInDeg(self):

        self.skySim.setStarRaDecInDeg(np.array([0]), np.array([1]),
                                      np.array([2]), np.array([3]))
        self.assertEqual(len(self.skySim.starId), 1)

    def testAddStarByFile(self):

        self._addStarByFile()
        self.assertEqual(len(self.skySim.starId), 8)
        self.assertEqual(self.skySim.ra[2], -1.176)
        self.assertEqual(self.skySim.decl[2], 1.196)
        self.assertEqual(self.skySim.mag[2], 17.0)

    def _addStarByFile(self):

        skyFile = os.path.join(getModulePath(), "tests", "testData", "sky",
                               "wfsStar.txt")
        self.skySim.addStarByFile(skyFile)

    def testExportSkyToFile(self):

        self._addStarByFile()
        outputFilePath = os.path.join(getModulePath(), "output",
                                      "testSkyOutput.txt")

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
        self.skySim.addStarByChipPos(sensorName, starId, xInpixelInCam,
                                     yInPixelInCam, starMag)

        # Test the result
        self.assertAlmostEqual(self.skySim.ra[0], 359.99971038)
        self.assertAlmostEqual(self.skySim.decl[0], 0.0001889)

    def _setObservationMetaData(self):

        ra = 0
        decl = 0
        rotSkyPos = 0
        mjd = 59580.0
        self.skySim.setObservationMetaData(ra, decl, rotSkyPos, mjd)

    def testGetCornOfChipOnSky(self):

        self._setObservationMetaData()

        sensorName = "R22_S11"
        cornerInRaDecList = self.skySim.getCornOfChipOnSky(sensorName)
        self.assertEqual(len(cornerInRaDecList), 4)

    def testAddStarByQueryDatabase(self):

        self._setDbInfo()

        # Query the database
        corner1 = [75.998622, -1]
        corner2 = [75.998622, -2]
        corner3 = [75.998985, -1]
        corner4 = [75.998985, -2]
        try:
            self.skySim.addStarByQueryDatabase(FilterType.U, corner1,
                                               corner2, corner3, corner4)
            # Check the adding of star
            self.assertEqual(len(self.skySim.starId), 3)
        except Exception as SystemExit:
            print("Please connect the UW database for the testing of query.")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
