import os
import numpy as np
import unittest

from lsst.sims.utils import ObservationMetaData
from lsst.obs.lsstSim import LsstSimMapper

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getModulePath


class TestSkySim(unittest.TestCase):
    """ Test the SkySim class."""

    def setUp(self):

        # Star file
        self.skyFile = os.path.join(getModulePath(), "tests", "testData",
                                    "sky", "wfsStar.txt")

        # Directory to the focal plane file
        self.testDataDir = os.path.join(getModulePath(), "tests", "testData",
                                        "testOpdFunc")

    def testFunc(self):

        # Instantiate the skySim object
        skySim = SkySim()

        skySim.setStarRaDecInDeg(np.array([0]), np.array([1]), np.array([2]), 
                                    np.array([3]))
        self.assertEqual(len(skySim.starId), 1)

        skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.assertEqual(len(skySim.starId), 2)

        # Check to add the second star with the same ID
        skySim.addStarByRaDecInDeg([1, 2], [2, 2], [3, 3], [4, 4])
        self.assertEqual(len(skySim.starId), 3)

        skySim.resetSky()
        self.assertEqual(len(skySim.starId), 0)

        skySim.addStarByFile(self.skyFile)
        self.assertEqual(len(skySim.starId), 8)

        outputFilePath = os.path.join(os.path.dirname(self.skyFile), "testSkyOutput.txt")
        if (not os.path.isfile(outputFilePath)):
            skySim.exportSkyToFile(outputFilePath)
            self.assertTrue(os.path.isfile(outputFilePath))
            os.remove(outputFilePath)
        else:
            print("Can't do the export file test because %s exists already." % outputFilePath)

    def testAddStarByChipPos(self):

        # Instantiate the skySim object
        skySim = SkySim()

        # Set the ObservationMetaData
        RA = 0
        Dec = 0
        cameraRotation = 0
        cameraMJD = 59580.0
        obs = ObservationMetaData(pointingRA=RA, pointingDec=Dec, rotSkyPos=cameraRotation, 
                                    mjd=cameraMJD)

        # Set the camera
        camera = LsstSimMapper().camera

        # Add the star
        sensorName = "R22_S11"
        starId = 0
        starMag = 17
        xInpixelInCam = 2000
        yInPixelInCam = 2036
        skySim.addStarByChipPos(camera, obs, sensorName, starId, xInpixelInCam, yInPixelInCam, 
                                starMag, self.testDataDir)

        # Test the result
        self.assertAlmostEqual(skySim.ra[0], 359.99971038)
        self.assertAlmostEqual(skySim.decl[0], 0.0001889)

        # Test to get the sensor box
        cornerInRaDecList = skySim.getCornOfChipOnSky(camera, obs, sensorName, self.testDataDir)
        self.assertEqual(len(cornerInRaDecList), 4)

    def testAddStarByQueryDatabase(self):

        # Instantiate the skySim object
        skySim = SkySim()

        # Filter type 
        aFilter = "u"

        # Remote database setting
        databaseHost = "localhost:51433"
        databaseUser = "LSST-2"
        databasePassword = "L$$TUser"
        databaseName = "LSSTCATSIM"

        # Config the database
        skySim.configDbInfo(databaseHost, databaseUser, databasePassword, databaseName)

        # Check the configuration
        self.assertEqual(skySim.dbInfo["host"], databaseHost)

        # Query the database
        corner1 = [75.998622, -1]
        corner2 = [75.998622, -2]
        corner3 = [75.998985, -1]
        corner4 = [75.998985, -2]
        try:
            skySim.addStarByQueryDatabase(aFilter, corner1, corner2, corner3, corner4)
            # Check the adding of star
            self.assertEqual(len(skySim.starId), 3)
        except Exception as SystemExit:
            print("Please connect the remote UW database for the testing of query.")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
