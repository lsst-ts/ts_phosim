import os, unittest
import numpy as np

from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.coordUtils.CameraUtils import raDecFromPixelCoords
from lsst.sims.utils import ObservationMetaData

from wep.SourceProcessor import SourceProcessor, expandDetectorName
from bsc.BrightStarDatabase import BrightStarDatabase

class SkySim(object):

    def __init__(self):
        """
        
        Initiate the SkySim object.
        """

        self.starId = np.array([])
        self.ra = np.array([])
        self.decl = np.array([])
        self.mag = np.array([])

        self.dbInfo = {"host":None, "user":None, "password":None, "dbName":None}

    def configDbInfo(self, host, user, password, dbName):
        """
        
        Configure the database information.
        
        Arguments:
            host {[str]} -- Database host name.
            user {[str]} -- Database user name.
            password {[str]} -- Database user password.
            dbName {[str]} -- Database name.
        """

        self.dbInfo["host"] = host
        self.dbInfo["user"] = user
        self.dbInfo["password"] = password
        self.dbInfo["dbName"] = dbName

    def setStarRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """
        
        Set the star information by (ra, dec) in degrees.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
            raInDeg {[list/ ndarray]} -- Star Ra in degree.
            declInDeg {[list/ ndarray]} -- Star Decl in degree.
            mag {[list/ ndarray]} -- Star magnitude.
        """

        # Clear the sky first
        self.resetSky()

        # Add the star
        self.addStarByRaDecInDeg(starId, raInDeg, declInDeg, mag)

    def __checkUniqStarId(self, starId):
        """
        
        Check the star ID is unique or not.
        
        Arguments:
            starId {[int]} -- Star Id.
        
        Returns:
            [bool] -- Is unique star Id or not.
        """

        isUnique = True
        if starId in self.starId:
            isUnique = False
            print("StarId=%d is not unique." % starId)

        return isUnique

    def addStarByRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """
        
        Add the star information by (ra, dec) in degrees.
        
        Arguments:
            starId {[int/ list/ ndarray]} -- Star ID.
            raInDeg {[float/ list/ ndarray]} -- Star Ra in degree.
            declInDeg {[float/ list/ ndarray]} -- Star Decl in degree.
            mag {[float/ list/ ndarray]} -- Star magnitude.
        """

        # Check the inputs are list or not, and change the type if necessary
        starId = self.__changeToList(starId)
        raInDeg = self.__changeToList(raInDeg)
        declInDeg = self.__changeToList(declInDeg)
        mag = self.__changeToList(mag)

        # Add the stars
        for ii in range(len(starId)):
            if (self.__checkUniqStarId(starId[ii])):
                self.starId = np.append(self.starId, int(starId[ii]))
                self.ra = np.append(self.ra, raInDeg[ii])
                self.decl = np.append(self.decl, declInDeg[ii])
                self.mag = np.append(self.mag, mag[ii])

    def __changeToList(self, variable):
        """
        
        Change the data type to list.
        
        Arguments:
            variable {[int/ float/ list/ ndarray]} -- Variable.
        
        Returns:
            [list] -- Variable as the list type.
        """

        try:
            len(variable)
        except Exception as TypeError:
            variable = [variable]

        return variable

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

        # Get the sky position in (ra, decl)
        raInDeg, declInDeg = self.__getSkyPosByChipPos(camera, obs, sensorName, xInpixelInCam, 
                                            yInPixelInCam, folderPath2FocalPlane, epoch=epoch, 
                                            includeDistortion=includeDistortion)

        # Add the star
        self.addStarByRaDecInDeg(starId, raInDeg, declInDeg, starMag)

    def getCornOfChipOnSky(self, camera, obs, sensorName, folderPath2FocalPlane, epoch=2000.0, 
                            includeDistortion=True):
        """
        
        Get the corner points of chip on sky position in (ra, dec). The direction is counter clockwise 
        from the origin.
        
        Arguments:
            camera {[Camera]} -- Camera object (e.g. LsstSimMapper().camera).
            obs {[ObservationMetaData]} -- Observation metadata object.
            sensorName {[str]} -- Abbreviated sensor name (e.g. "R22_S11").
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
        
        Returns:
            [list] -- List of corner points in tuple. The unit is in (ra, dec).
        """

        # Get the sensor dimension
        sourProc = SourceProcessor()
        sourProc.config(folderPath2FocalPlane=folderPath2FocalPlane)
        dimX, dimY = sourProc.sensorDimList[sensorName]

        # Corner in Pixel
        xInPixel = [0, dimX, dimX, 0]
        yInPixel = [0, 0, dimY, dimY]

        # Corner in (Ra, Dec)
        cornerInRaDecList = []
        for ii in range(4):
            raInDeg, declInDeg = self.__getSkyPosByChipPos(camera, obs, sensorName, xInPixel[ii], 
                                                yInPixel[ii], folderPath2FocalPlane, epoch=epoch, 
                                                includeDistortion=includeDistortion)
            cornerInRaDecList.append((raInDeg, declInDeg))

        return cornerInRaDecList

    def __getSkyPosByChipPos(self, camera, obs, sensorName, xInpixelInCam, yInPixelInCam, 
                                folderPath2FocalPlane, epoch=2000.0, includeDistortion=True):
        """
        
        Get the sky position in (ra, dec) based on the chip pixel positions.
        
        Arguments:
            camera {[Camera]} -- Camera object (e.g. LsstSimMapper().camera).
            obs {[ObservationMetaData]} -- Observation metadata object.
            sensorName {[str]} -- Abbreviated sensor name (e.g. "R22_S11").
            xInpixelInCam {[float]} -- Pixel position x on camera coordinate.
            yInPixelInCam {[float]} -- Pixel position y on camera coordinate.
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
        
        Returns:
            [float] -- Ra in degree.
            [float] -- Decl in degree.
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

        return raInDeg, declInDeg

    def addStarByQueryDatabase(self, aFilter, corner1, corner2, corner3, corner4, 
                                tableName="bright_stars"):
        """
        
        Add the star by querying the database of bright star catalog in UW.
        
        Arguments:
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
            corner1 {[float]} -- The first corner of the sensor defined as (RA, Decl).
            corner2 {[float]} -- The second corner of the sensor defined as (RA, Decl).
            corner3 {[float]} -- The third corner of the sensor defined as (RA, Decl).
            corner4 {[float]} -- The fourth corner of the sensor defined as (RA, Decl).
        
        Keyword Arguments:
            tableName {str} -- Table name in database. (default: {"bright_stars"})
        """

        # Instantiate the BrightStarDatabase
        bsc = BrightStarDatabase()

        # Connect to the database
        bsc.connect(self.dbInfo["host"], self.dbInfo["user"], self.dbInfo["password"], 
                    self.dbInfo["dbName"])

        # Query the database        
        stars = bsc.query(tableName, aFilter, corner1, corner2, corner3, corner4)

        # Disconnect from the database
        bsc.disconnect()

        # Add the star information
        starId = np.array(stars.SimobjID).astype("int")
        mag = getattr(stars, "LSSTMag"+aFilter.upper())
        self.addStarByRaDecInDeg(starId, stars.RA, stars.Decl, mag)

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

    # Do the unit test
    unittest.main()
