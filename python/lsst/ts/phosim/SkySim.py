import numpy as np

from lsst.obs.lsstSim import LsstSimMapper
from lsst.sims.utils import ObservationMetaData
from lsst.sims.coordUtils.CameraUtils import raDecFromPixelCoords

from lsst.ts.wep.SourceProcessor import SourceProcessor, expandDetectorName
from lsst.ts.wep.bsc.BrightStarDatabase import BrightStarDatabase


class SkySim(object):

    def __init__(self):
        """Initialization of sky simulator class."""

        self.starId = np.array([], dtype=int)
        self.ra = np.array([])
        self.decl = np.array([])
        self.mag = np.array([])

        self.dbInfo = {"host": None, "user": None, "password": None,
                       "dbName": None}

        self._camera = LsstSimMapper().camera
        self._obs = ObservationMetaData()
        self._sourProc = SourceProcessor()

    def setCamera(self, camera):
        """Set the camera object.

        Parameters
        ----------
        camera : Camera
            A collection of Detectors that also supports coordinate
            transformation
        """

        self._camera = camera

    def setObservationMetaData(self, ra, decl, rotSkyPos, mjd):
        """Set the observation meta data.

        Parameters
        ----------
        ra : float
            Pointing ra in degree.
        decl : float
            Pointing decl in degree.
        rotSkyPos : float
            The orientation of the telescope in degrees.
        mjd : float
            Camera MJD.
        """

        self._obs = ObservationMetaData(pointingRA=ra, pointingDec=decl,
                                        rotSkyPos=rotSkyPos,
                                        mjd=mjd)

    def setFolderPath2FocalPlane(self, folderPath2FocalPlane):
        """Set the folder path to focal plane data.

        Parameters
        ----------
        folderPath2FocalPlane : str
            Folder path to focal plane data.
        """

        self._sourProc.config(folderPath2FocalPlane=folderPath2FocalPlane)

    def configDbInfo(self, host, user, password, dbName):
        """Configure the database information.

        Parameters
        ----------
        host : str
            Database host name.
        user : str
            Database user name.
        password : str
            Database user password.
        dbName : str
            Database name.
        """

        self.dbInfo["host"] = host
        self.dbInfo["user"] = user
        self.dbInfo["password"] = password
        self.dbInfo["dbName"] = dbName

    def addStarByRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """Add the star information by (ra, dec) in degrees.

        Parameters
        ----------
        starId : int, list[int], or numpy.ndarray[int]
            Star Id.
        raInDeg : float, list, or numpy.ndarray
            Star ra in degree.
        declInDeg : float, list, or numpy.ndarray
            Star decl in degree.
        mag : float, list, or numpy.ndarray
            Star magnitude.
        """

        # Check the inputs are list or not, and change the type if necessary
        starIdList = self._changeToListIfNecessary(starId)
        raInDegList = self._changeToListIfNecessary(raInDeg)
        declInDegList = self._changeToListIfNecessary(declInDeg)
        magList = self._changeToListIfNecessary(mag)

        # Add the stars
        for ii in range(len(starIdList)):
            intStarId = int(starIdList[ii])
            if (self._isUniqStarId(intStarId)):
                self.starId = np.append(self.starId, intStarId)
                self.ra = np.append(self.ra, raInDegList[ii])
                self.decl = np.append(self.decl, declInDegList[ii])
                self.mag = np.append(self.mag, magList[ii])

    def _changeToListIfNecessary(self, variable):
        """Change the data type to list.

        Parameters
        ----------
        variable : int, float, list, or numpy.ndarray
            Variable.

        Returns
        -------
        list
            Variable as the list.
        """

        if isinstance(variable, (int, float)):
            return [variable]
        else:
            return variable

    def _isUniqStarId(self, starId):
        """Check the star ID is unique or not.

        Parameters
        ----------
        starId : int
            Star Id.

        Returns
        -------
        bool
            True if the unique Id.
        """

        if starId in self.starId:
            isUnique = False
            print("StarId=%d is not unique." % starId)
        else:
            isUnique = True

        return isUnique

    def resetSky(self):
        """Reset the sky information and delete all existed stars."""

        self.__init__()

    def setStarRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """Set the star information by (ra, dec) in degrees.

        Parameters
        ----------
        starId : int, list[int], or numpy.ndarray[int]
            Star Id.
        raInDeg : float, list, or numpy.ndarray
            Star ra in degree.
        declInDeg : float, list, or numpy.ndarray
            Star decl in degree.
        mag : float, list, or numpy.ndarray
            Star magnitude.
        """

        self.resetSky()
        self.addStarByRaDecInDeg(starId, raInDeg, declInDeg, mag)

    def addStarByFile(self, readFilePath, skiprows=0):
        """Add the star data by reading the file.

        Parameters
        ----------
        readFilePath : str
            Star data file path.
        skiprows : int, optional
            Skip the first "skiprows" lines. (the default is 0.)
        """

        data = np.loadtxt(readFilePath, skiprows=skiprows)
        for star in data:
            self.addStarByRaDecInDeg(star[0], star[1], star[2], star[3])

    def exportSkyToFile(self, outputFilePath):
        """Export the star information into the file.

        Parameters
        ----------
        outputFilePath : str
            Output file path.
        """

        # Add the header (star ID, ra, decl, magnitude)
        content = "# Id\t Ra\t\t Decl\t\t Mag\n"

        # Add the star information
        for ii in range(len(self.starId)):
            content += "%d\t %3.6f\t %3.6f\t %3.6f\n" % (
                        self.starId[ii], self.ra[ii], self.decl[ii],
                        self.mag[ii])

        # Write into file
        fid = open(outputFilePath, "w")
        fid.write(content)
        fid.close()

    def addStarByChipPos(self, sensorName, starId, xInpixelInCam,
                         yInPixelInCam, starMag, epoch=2000.0,
                         includeDistortion=True):
        """Add the star based on the chip position.

        Parameters
        ----------
        sensorName : str
            Abbreviated sensor name (e.g. "R22_S11").
        starId : int
            Star Id.
        xInpixelInCam : float
            Pixel position x on camera coordinate.
        yInPixelInCam : float
            Pixel position y on camera coordinate.
        starMag : float
            Star magnitude.
        epoch : float, optional
            Epoch is the mean epoch in years of the celestial coordinate
            system. (the default is 2000.0.)
        includeDistortion : bool, optional
            If True (default), then this method will expect the true pixel
            coordinates with optical distortion included.  If False, this
            method will expect TAN_PIXEL coordinates, which are the pixel
            coordinates with estimated optical distortion removed.  See
            the documentation in afw.cameraGeom for more details. (the
            default is True.)
        """

        # Get the sky position in (ra, decl)
        raInDeg, declInDeg = self._getSkyPosByChipPos(
                                      sensorName, xInpixelInCam,
                                      yInPixelInCam, epoch=epoch,
                                      includeDistortion=includeDistortion)

        # Add the star
        self.addStarByRaDecInDeg(starId, raInDeg, declInDeg, starMag)

    def getCornOfChipOnSky(self, sensorName, epoch=2000.0,
                           includeDistortion=True):
        """Get the corner points of chip on sky position in (ra, dec).

        The direction is counter clockwise from the origin.

        Parameters
        ----------
        sensorName : str
            Abbreviated sensor name (e.g. "R22_S11").
        epoch : float, optional
            Epoch is the mean epoch in years of the celestial coordinate
            system. (the default is 2000.0.)
        includeDistortion : bool, optional
            If True (default), then this method will expect the true pixel
            coordinates with optical distortion included.  If False, this
            method will expect TAN_PIXEL coordinates, which are the pixel
            coordinates with estimated optical distortion removed.  See
            the documentation in afw.cameraGeom for more details. (the
            default is True.)

        Returns
        -------
        list
            List of corner points in tuple. The unit is in (ra, dec).
        """

        # Get the sensor dimension
        dimX, dimY = self._sourProc.sensorDimList[sensorName]

        # Corner in Pixel
        xInPixel = [0, dimX, dimX, 0]
        yInPixel = [0, 0, dimY, dimY]

        # Corner in (Ra, Dec)
        cornerInRaDecList = []
        for ii in range(4):
            raInDeg, declInDeg = self._getSkyPosByChipPos(
                                        sensorName, xInPixel[ii],
                                        yInPixel[ii], epoch=epoch,
                                        includeDistortion=includeDistortion)
            cornerInRaDecList.append((raInDeg, declInDeg))

        return cornerInRaDecList

    def _getSkyPosByChipPos(self, sensorName, xInpixelInCam, yInPixelInCam,
                            epoch=2000.0, includeDistortion=True):
        """Get the sky position in (ra, dec) based on the chip pixel positions.

        Parameters
        ----------
        sensorName : str
            Abbreviated sensor name (e.g. "R22_S11").
        xInpixelInCam : float
            Pixel position x on camera coordinate.
        yInPixelInCam : float
            Pixel position y on camera coordinate.
        epoch : float, optional
            Epoch is the mean epoch in years of the celestial coordinate
            system. (the default is 2000.0.)
        includeDistortion : bool, optional
            If True (default), then this method will expect the true pixel
            coordinates with optical distortion included.  If False, this
            method will expect TAN_PIXEL coordinates, which are the pixel
            coordinates with estimated optical distortion removed.  See
            the documentation in afw.cameraGeom for more details. (the
            default is True.)

        Returns
        -------
        float
            Ra in degree.
        float
            Decl in degree.
        """

        # Get the pixel positions in DM team
        self._sourProc.config(sensorName=sensorName)
        pixelDmX, pixelDmY = self._sourProc.camXY2DmXY(xInpixelInCam,
                                                       yInPixelInCam)

        # Expend the sensor name
        expendedSensorName = expandDetectorName(sensorName)

        # Get the sky position in (ra, decl)
        raInDeg, declInDeg = raDecFromPixelCoords(
                                pixelDmX, pixelDmY, expendedSensorName,
                                camera=self._camera, obs_metadata=self._obs,
                                epoch=epoch,
                                includeDistortion=includeDistortion)

        return raInDeg, declInDeg

    def addStarByQueryDatabase(self, filterType, corner1, corner2, corner3,
                               corner4, tableName="bright_stars"):
        """Add the star by querying the database of bright star catalog in UW.

        Parameters
        ----------
        filterType : FilterType
            Active filter type.
        corner1 : float
            The first corner of the sensor defined as (RA, Decl).
        corner2 : float
            The second corner of the sensor defined as (RA, Decl).
        corner3 : float
            The third corner of the sensor defined as (RA, Decl).
        corner4 : float
            The fourth corner of the sensor defined as (RA, Decl).
        tableName : str, optional
            Table name in database. (the default is "bright_stars".)
        """

        # Instantiate the BrightStarDatabase
        bsc = BrightStarDatabase()

        # Connect to the database
        bsc.connect(self.dbInfo["host"], self.dbInfo["user"],
                    self.dbInfo["password"], self.dbInfo["dbName"])

        # Query the database
        aFilter = filterType.name.lower()
        stars = bsc.query(tableName, aFilter, corner1, corner2, corner3,
                          corner4)

        # Disconnect from the database
        bsc.disconnect()

        # Add the star information
        starId = np.array(stars.SimobjID).astype("int")
        mag = getattr(stars, "LSSTMag"+aFilter.upper())
        self.addStarByRaDecInDeg(starId, stars.RA, stars.Decl, mag)


if __name__ == "__main__":
    pass