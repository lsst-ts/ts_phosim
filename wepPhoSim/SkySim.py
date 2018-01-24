import os, unittest
import numpy as np

from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils.CameraUtils import raDecFromPixelCoords
from lsst.sims.utils import ObservationMetaData

from wep.SourceProcessor import SourceProcessor, expandDetectorName

class SkySim(object):

    def __init__(self):
        """
        
        Initiate the SkySim object.
        """

        self.starId = np.array([])
        self.ra = np.array([])
        self.decl = np.array([])
        self.mag = np.array([])

    def setStarRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """
        
        Set the star information by (ra, dec) in degrees.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
            raInDeg {[list/ ndarray]} -- Star Ra in degree.
            declInDeg {[list/ ndarray]} -- Star Decl in degree.
            mag {[list/ ndarray]} -- Star magnitude.
        """

        self.__setUniqStarId(starId)
        self.ra = raInDeg
        self.decl = declInDeg
        self.mag = mag

    def __setUniqStarId(self, starId):
        """
        
        Set the star unique IDs.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
        
        Raises:
            ValueError -- Star IDs are not unique.
        """

        # Collect all star id
        allStarId = np.append(self.starId, starId).astype("int")

        # Check all IDs are unique all not
        if (len(allStarId) != len(np.unique(allStarId))):
            raise ValueError("Star IDs are not unique.")

        # Set the star ID
        self.starId = allStarId

    def addStarByRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """
        
        Add the star information by (ra, dec) in degrees.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
            raInDeg {[list/ ndarray]} -- Star Ra in degree.
            declInDeg {[list/ ndarray]} -- Star Decl in degree.
            mag {[list/ ndarray]} -- Star magnitude.
        """

        self.__setUniqStarId(starId)
        self.ra = np.append(self.ra, raInDeg)
        self.decl = np.append(self.decl, declInDeg)
        self.mag = np.append(self.mag, mag)

    def resetSky(self):
        """
        
        Reset the sky information and delete all existed stars.
        """
        
        self.__init__()

    def exportSkyToFile(self, outputFilePath):
        """
        
        Export the star information into the file.
        
        Arguments:
            outputFilePath {[str]} -- Output file path.
        """
        
        # Add the header (star ID, ra, decl, magnitude)
        content = "Id\t Ra\t\t Decl\t\t Mag\n"

        # Add the star information
        for ii in range(len(self.starId)):
            content += "%d\t %3.6f\t %3.6f\t %3.6f\n" % (self.starId[ii], self.ra[ii], 
                                                         self.decl[ii], self.mag[ii])
        # Write into file
        fid = open(outputFilePath, "w")
        fid.write(content)
        fid.close()

    def addStarByFile(self, readFilePath, skiprows=1):
        """
        
        Add the star data by reading the file.
        
        Arguments:
            readFilePath {[star]} -- Star data file path.
        
        Keyword Arguments:
            skiprows {int} -- Skip the first "skiprows" lines. (default: {1})
        """
        
        # Read the data
        data = np.loadtxt(readFilePath, skiprows=skiprows)

        # Add the stars
        for star in data:
            self.addStarByRaDecInDeg(star[0], star[1], star[2], star[3])

    def addStarByChipPos(self, camera, obs, sensorName, starId, xInpixelInCam, yInPixelInCam, 
                            starMag, folderPath2FocalPlane, epoch=2000.0, includeDistortion=True):
        """
        
        Add the star based on the chip position.
        
        Arguments:
            camera {[Camera]} -- Camera object (e.g. LsstSimMapper().camera).
            obs {[ObservationMetaData]} -- Observation metadata object.
            sensorName {[str]} -- Abbreviated sensor name (e.g. "R22_S11").
            starId {[int]} -- Star ID.
            xInpixelInCam {[float]} -- Pixel position x on camera coordinate.
            yInPixelInCam {[float]} -- Pixel position y on camera coordinate.
            starMag {[float]} -- Star magnitude.
            folderPath2FocalPlane {[str]} -- Path to the directory of focal plane data 
                                            ("focalplanelayout.txt").
        
        Keyword Arguments:
            epoch {float} -- Epoch is the mean epoch in years of the celestial coordinate system. 
                            (default: {2000.0})
            includeDistortion {bool} -- If True (default), then this method will expect the true pixel 
                                        coordinates with optical distortion included.  If False, this 
                                        method will expect TAN_PIXEL coordinates, which are the pixel 
                                        coordinates with estimated optical distortion removed.  See 
                                        the documentation in afw.cameraGeom for more details. 
                                        (default: {True})
        """

        # Get the pixel positions in DM team
        sourProc = SourceProcessor()
        sourProc.config(sensorName=sensorName, folderPath2FocalPlane=folderPath2FocalPlane)
        pixelDmX, pixelDmY = sourProc.camXY2DmXY(xInpixelInCam, yInPixelInCam)

        # Expend the sensor name
        expendedSensorName = expandDetectorName(sensorName)

        # Get the sky position in (ra, decl)
        raInDeg, declInDeg = raDecFromPixelCoords(pixelDmX, pixelDmY, expendedSensorName, 
                                                camera=camera, obs_metadata=obs, epoch=epoch, 
                                                includeDistortion=includeDistortion)

        # Add the star
        self.addStarByRaDecInDeg(starId, raInDeg, declInDeg, starMag)

class SkySimTest(unittest.TestCase):
    
    """
    Test functions in SkySim.
    """

    def setUp(self):

        # Star file
        self.skyFile = os.path.join("..", "data", "sky", "wfsStar.txt")

        # Directory to the focal plane file
        self.testDataDir = os.path.join("..", "testData", "testOpdFunc")

    def testFunc(self):

        # Instantiate the skySim object
        skySim = SkySim()

        skySim.setStarRaDecInDeg(np.array([0]), np.array([1]), np.array([2]), 
                                    np.array([3]))
        self.assertEqual(len(skySim.starId), 1)

        skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.assertEqual(len(skySim.starId), 2)

        try:
            skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        except Exception as ValueError:
            pass

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

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
