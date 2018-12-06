import numpy as np
from astropy.io import fits

from lsst.ts.wep.cwfs.Tool import ZernikeAnnularFit, ZernikeEval
from lsst.ts.wep.SourceProcessor import SourceProcessor
from lsst.ts.phosim.MetroTool import calc_pssn, psf2eAtmW


class OpdMetrology(object):

    def __init__(self):
        """Initialization of OPD metrology class.

        OPD: Optical path difference.
        """

        self.wt = np.array([])
        self.fieldX = np.array([])
        self.fieldY = np.array([])

    def setWeightingRatio(self, wt):
        """Set the weighting ratio used in Gaussian quadrature.

        Parameters
        ----------
        wt : list or numpy.ndarray
            Weighting ratio.
        """

        wtArray = np.array(wt, dtype=float)
        self._setNormalizedWeightingRatio(wtArray)

    def _setNormalizedWeightingRatio(self, wtArray):
        """Normalize the weighting ratios.

        Parameters
        ----------
        wtArray : numpy.ndarray
            Weighting ratio.

        Raises
        ------
        ValueError
            All weighting ratios should be >=0.
        """

        if np.all(wtArray >= 0):
            self.wt = wtArray/np.sum(wtArray)
        else:
            raise ValueError("All weighting ratios should be >=0.")

    def setFieldXYinDeg(self, fieldXInDegree, fieldYInDegree):
        """Set the field X, Y in degree.

        Parameters
        ----------
        fieldXInDegree : list or numpy.ndarray
            Field X in degree.
        fieldYInDegree : list or numpy.ndarray
            Field Y in degree.
        """

        self.fieldX = np.array(fieldXInDegree, dtype=float)
        self.fieldY = np.array(fieldYInDegree, dtype=float)

    def addFieldXYbyDeg(self, fieldXInDegree, fieldYInDegree):
        """Add the new field X, Y in degree.

        Parameters
        ----------
        fieldXInDegree : float, list, or numpy.ndarray
            New field X in degree.
        fieldYInDegree : float, list, or numpy.ndarray
            New field Y in degree.
        """

        self.fieldX = np.append(self.fieldX, fieldXInDegree)
        self.fieldY = np.append(self.fieldY, fieldYInDegree)

    def setDefaultLsstGQ(self):
        """Set the default LSST GQ field X, Y and weighting ratio.

        GQ: Gaussian quadrature
        """

        # The distance of point xi (used in Gaussian quadrature plane) to the
        # origin
        # This value is in [-1.75, 1.75]
        armLen = [0.379, 0.841, 1.237, 1.535, 1.708]

        # Weighting of point xi (used in Gaussian quadrature plane) for each
        # ring
        armW = [0.2369, 0.4786, 0.5689, 0.4786, 0.2369]

        # Number of points on each ring
        nArm = 6

        # Get the weighting for all field points (31 for lsst camera)
        # Consider the first element is center (0)
        wt = np.concatenate([np.zeros(1), np.kron(armW, np.ones(nArm))])
        self.setWeightingRatio(wt)

        # Generate the fields point x, y coordinates
        pointAngle = np.arange(nArm) * (2*np.pi)/nArm
        fieldX = np.concatenate([np.zeros(1),
                                 np.kron(armLen, np.cos(pointAngle))])
        fieldY = np.concatenate([np.zeros(1),
                                 np.kron(armLen, np.sin(pointAngle))])
        self.setFieldXYinDeg(fieldX, fieldY)

    def getDefaultLsstWfsGQ(self):
        """Get the default field X, Y of LSST WFS on GQ.

        WFS: Wavefront sensor.
        GQ: Gaussian quadrature

        Returns
        -------
        list
            Field X in degree.
        list
            Field Y in degree.
        """

        # Field x, y for 4 WFS
        fieldWFSx = [1.176, -1.176, -1.176, 1.176]
        fieldWFSy = [1.176, 1.176, -1.176, -1.176]

        return fieldWFSx, fieldWFSy

    def setDefaultComcamGQ(self):
        """Set the default ComCam GQ field X, Y and weighting ratio.

        GQ: Gaussian quadrature
        """

        # ComCam is the cetral raft of LSST cam, which is composed of 3 x 3
        # CCDs.
        nRow = 3
        nCol = 3

        # Number of field points
        nField = nRow*nCol

        # Get the weighting for all field points (9 for comcam)
        wt = np.ones(nField)
        self.setWeightingRatio(wt)

        # Distance to raft center in degree along x/y direction and the
        # related relative position
        sensorD = 0.2347
        coorComcam = sensorD * np.array([-1, 0, 1])

        # Generate the fields point x, y coordinates
        fieldX = np.kron(coorComcam, np.ones(nRow))
        fieldY = np.kron(np.ones(nCol), coorComcam)
        self.setFieldXYinDeg(fieldX, fieldY)

    def getZkFromOpd(self, opdFitsFile=None, opdMap=None, znTerms=22,
                     obscuration=0.61):
        """Get the wavefront error of OPD in the basis of annular Zernike
        polynomials.

        OPD: Optical path difference.

        Parameters
        ----------
        opdFitsFile : str, optional
            OPD FITS file. (the default is None)
        opdMap : numpy.ndarray, optional
            OPD map data. (the default is None, which [default_description])
        znTerms : int, optional
            Number of terms of annular Zk (z1-z22 by default). (the default
            is 22.)
        obscuration : float, optional
            Obscuration of annular Zernike polynomial. (the default is 0.61.)

        Returns
        -------
        numpy.ndarray
            Annular Zernike polynomials. For PhoSim OPD, the unit is um.
        numpy.ndarray
            OPD map.
        numpy.ndarray
            Meshgrid x in OPD map.
        numpy.ndarray
            Meshgrid y in OPD map.

        Raises
        ------
        ValueError
            The x, y dimensions of OPD are different.
        """

        # Get the OPD data (PhoSim OPD unit: um)
        if (opdFitsFile is not None):
            opd = fits.getdata(opdFitsFile)
        elif (opdMap is not None):
            opd = opdMap.copy()

        # Check the x, y dimensions of OPD are the same
        if (np.unique(opd.shape).size != 1):
            raise ValueError("The x, y dimensions of OPD are different.")

        # x-, y-coordinate in the OPD image
        opdSize = opd.shape[0]
        opdGrid1d = np.linspace(-1, 1, opdSize)
        opdx, opdy = np.meshgrid(opdGrid1d, opdGrid1d)

        # Fit the OPD map with Zk and write into the file
        idx = (opd != 0)
        zk = ZernikeAnnularFit(opd[idx], opdx[idx], opdy[idx], znTerms,
                               obscuration)

        return zk, opd, opdx, opdy

    def rmPTTfromOPD(self, opdFitsFile=None, opdMap=None):
        """Remove the afftection of piston (z1), x-tilt (z2), and y-tilt (z3)
        from the OPD map.

        OPD: Optical path difference.

        Parameters
        ----------
        opdFitsFile : str, optional
            OPD FITS file. (the default is None.)
        opdMap : numpy.ndarray, optional
            OPD map data. (the default is None.)

        Returns
        -------
        numpy.ndarray
            OPD map after removing the affection of z1-z3.
        numpy.ndarray
            Meshgrid x in OPD map.
        numpy.ndarray
            Meshgrid y in OPD map.
        """

        # Do the spherical Zernike fitting for the OPD map
        # Only fit the first three terms (z1-z3): piston, x-tilt, y-tilt
        zk, opd, opdx, opdy = self.getZkFromOpd(
                                    opdFitsFile=opdFitsFile, opdMap=opdMap,
                                    znTerms=3, obscuration=0)

        # Find the index that the value of OPD is not 0
        idx = (opd != 0)

        # Remove the PTT
        opd[idx] -= ZernikeEval(zk, opdx[idx], opdy[idx])

        return opd, opdx, opdy

    def addFieldXYbyCamPos(self, sensorName, xInpixel, yInPixel,
                           folderPath2FocalPlane):
        """Add the new field X, Y in degree by the camera pixel positions.

        Parameters
        ----------
        sensorName : str
            Canera sensor name (e.g. "R22_S11", "R40_S02_C0").
        xInpixel : float
            Pixel x on camera coordinate.
        yInPixel : float
            Pixel y on camera coordinate.
        folderPath2FocalPlane : str
            Path to the directory of focal plane data ("focalplanelayout.txt").
        """

        # Get the focal plane data and set the sensor name
        sourProc = SourceProcessor()
        sourProc.config(sensorName=sensorName,
                        folderPath2FocalPlane=folderPath2FocalPlane)

        # Do the coordinate transformation
        fieldXInDegree, fieldYInDegree = sourProc.camXYtoFieldXY(
                                                    xInpixel, yInPixel)

        # Add to listed field x, y
        self.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

    def calcPSSN(self, wavelengthInUm, opdFitsFile=None, opdMap=None, zen=0,
                 debugLevel=0):
        """ Calculate the PSSN based on OPD map.

        PSSN: Normalized point source sensitivity.
        OPD: Optical path difference.

        Parameters
        ----------
        wavelengthInUm : float
            Wavelength in microns.
        opdFitsFile : str, optional
            OPD FITS file. (the default is None.)
        opdMap : numpy.ndarray, optional
            OPD map data. (the default is None.)
        zen : float, optional
            elescope zenith angle in degree. (the default is 0.)
        debugLevel : int, optional
            Debug level. The higher value gives more information. (the default
            is 0.)

        Returns
        -------
        float
            Calculated PSSN.
        """

        # Before calc_pssn,
        # (1) Remove PTT (piston, x-tilt, y-tilt),
        # (2) Make sure outside of pupil are all zeros
        opdRmPTT = self.rmPTTfromOPD(opdFitsFile=opdFitsFile, opdMap=opdMap)[0]

        # Calculate the normalized point source sensitivity (PSSN)
        pssn = calc_pssn(opdRmPTT, wavelengthInUm, zen=zen,
                         debugLevel=debugLevel)

        return pssn

    def calcFWHMeff(self, pssn):
        """Calculate the effective FWHM.

        FWHM: Full width at half maximum.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn : float
            PSSN value.

        Returns
        -------
        float
            Effective FWHM.
        """

        # FWHMeff_sys = FWHMeff_atm * sqrt(1/PSSN - 1).
        # FWHMeff_atm = 0.6 arcsec.
        # Another correction factor (eta = 1.086) is used to account for the
        # difference between the
        # simple RSS and the more proper convolution.
        # Follow page 7 (section 7.2 FWHM) in document-17242 for more
        # information.
        eta = 1.086
        FWHMatm = 0.6
        FWHMeff = eta*FWHMatm*np.sqrt(1/pssn - 1)

        return FWHMeff

    def calcDm5(self, pssn):
        """Calculate the loss of limiting depth (dm5).

        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssn : float
            PSSN value.

        Returns
        -------
        float
            dm5.
        """

        # Calculate dm5 (the loss of limiting depth)
        # Check eq. (4.1) in page 4 in document-17242 for more information.
        dm5 = -1.25 * np.log10(pssn)

        return dm5

    def calcEllip(self, wavelengthInUm, opdFitsFile=None, opdMap=None, zen=0,
                  debugLevel=0):
        """Calculate the ellipticity.

        Parameters
        ----------
        wavelengthInUm : float
            Wavelength in microns.
        opdFitsFile : str, optional
            OPD FITS file. (the default is None.)
        opdMap : numpy.ndarray, optional
            OPD map data. (the default is None.)
        zen : float, optional
            elescope zenith angle in degree. (the default is 0.)
        debugLevel : int, optional
            Debug level. The higher value gives more information. (the default
            is 0.)

        Returns
        -------
        float
            Ellipticity.
        """

        # Remove the affection of piston (z1), x-tilt (z2), and y-tilt (z3)
        # from OPD map.
        opdRmPTT = self.rmPTTfromOPD(opdFitsFile=opdFitsFile, opdMap=opdMap)[0]

        # Calculate the ellipticity
        elli = psf2eAtmW(opdRmPTT, wavelengthInUm, zen=zen,
                         debugLevel=debugLevel)[0]

        return elli

    def calcGQvalue(self, valueList):
        """Calculate the GQ value.

        GQ: Gaussian quadrature

        Parameters
        ----------
        valueList : list or numpy.ndarray
            List of value (PSSN, effective FWHM, dm5, ellipticity).

        Returns
        -------
        float
            GQ value.

        Raises
        ------
        ValueError
            Length of wt ratio != length of value list.
        """

        # Check the lengths of weighting ratio and value list are the same
        if (len(self.wt) != len(valueList)):
            raise ValueError("Length of wt ratio != length of value list.")

        # Calculate the effective value on Gaussain quardure plane
        valueArray = np.array(valueList, dtype=float)
        GQvalue = np.sum(self.wt*valueArray)

        return GQvalue


if __name__ == "__main__":
    pass
