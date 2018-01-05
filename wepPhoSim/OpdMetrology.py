import numpy as np
from astropy.io import fits

from cwfs.Tool import ZernikeAnnularFit, ZernikeAnnularEval, ZernikeEval
from wep.SourceProcessor import SourceProcessor
from wepPhoSim.MetroTool import calc_pssn, psf2eAtmW

class OpdMetrology(object):

    def __init__(self):
        """
        
        Initiate the OpdMetrology object.
        """
        
        self.wt = np.array([])
        self.fieldX = np.array([])
        self.fieldY = np.array([])

    def setWeightingRatio(self, wt):
        """
        
        Set the weighting ratio used in Gaussian quadrature.
        
        Arguments:
            wt {[ndarray]} -- Weighting ratio.
        """

        self.wt = wt

    def setFieldXYinDeg(self, fieldXInDegree, fieldYInDegree):
        """
        
        Set the field X, Y in degree.
        
        Arguments:
            fieldXInDegree {[ndarray]} -- Field X in degree.
            fieldYInDegree {[ndarray]} -- Field Y in degree.
        """

        self.fieldX = fieldXInDegree
        self.fieldY = fieldYInDegree

    def addWeightingRatio(self, newWt):
        """
        
        Add the new weighting ratio values used in Gaussian quadrature..
        
        Arguments:
            newWt {[float/ list/ ndarray]} -- New weighting ratio.
        """
        
        self.wt = np.append(self.wt, newWt)

    def normalizeWeightingRatio(self):
        """
        
        Normalize the weighting ratios.
        
        Raises:
            ValueError -- Not all weighting ratios >=0.
        """

        # Check all values >=0
        if np.any(self.wt < 0):
            raise ValueError("All weighting ratios should be >=0.")

        # Do the normalization
        self.wt = self.wt/np.sum(self.wt)

    def addFieldXYbyDeg(self, fieldXInDegree, fieldYInDegree):
        """
        
        Add the new field X, Y in degree.
        
        Arguments:
            fieldXInDegree {[float/ list/ ndarray]} -- New field X in degree.
            fieldYInDegree {[float/ list/ ndarray]} -- New field Y in degree.
        """

        self.fieldX = np.append(self.fieldX, fieldXInDegree)
        self.fieldY = np.append(self.fieldY, fieldYInDegree)

    def setDefaultLsstGQ(self):
        """
        
        Set the default LSST Gaussian quadrature (GQ) field X, Y and weighting ratio.
        """

        # The distance of point xi (used in Gaussian quadrature plane) to the origin
        # This value is in [-1.75, 1.75]
        armLen = [0.379, 0.841, 1.237, 1.535, 1.708]

        # Weighting of point xi (used in Gaussian quadrature plane) for each ring
        armW = [0.2369, 0.4786, 0.5689, 0.4786, 0.2369]

        # Number of points on each ring
        nArm = 6

        # Get the weighting for all field points (31 for lsst camera)
        # Consider the first element is center (0)
        wt = np.concatenate([np.zeros(1), np.kron(armW, np.ones(nArm))])
        self.setWeightingRatio(wt)
        self.normalizeWeightingRatio()

        # Generate the fields point x, y coordinates
        pointAngle = np.arange(nArm) * (2*np.pi)/nArm
        fieldX = np.concatenate([np.zeros(1), np.kron(armLen, np.cos(pointAngle))])
        fieldY = np.concatenate([np.zeros(1), np.kron(armLen, np.sin(pointAngle))])
        self.setFieldXYinDeg(fieldX, fieldY)

    def getDefaultLsstWfsGQ(self):
        """
        
        Get the default field X, Y of LSST wavefront sensor (WFS) on Gaussian quadrature (GQ).
        
        Returns:
            [list] -- Field X in degree.
            [list] -- Field Y in degree.
        """

        # Field x, y for 4 WFS
        fieldWFSx = [1.176, -1.176, -1.176, 1.176]
        fieldWFSy = [1.176, 1.176, -1.176, -1.176]

        return fieldWFSx, fieldWFSy

    def setDefaultComcamGQ(self):
        """
        
        Set the default ComCam Gaussian quadrature (GQ) field X, Y and weighting ratio.
        """

        # ComCam is the cetral raft of LSST cam, which is composed of 3 x 3 CCDs.
        nRow = 3
        nCol = 3
        
        # Number of field points
        nField = nRow*nCol

        # Number of field points with the consideration of wavefront sensor
        nFieldWfs = nField

        # Get the weighting for all field points (9 for comcam)
        wt = np.ones(nField)
        self.setWeightingRatio(wt)
        self.normalizeWeightingRatio()

        # Distance to raft center in degree along x/y direction and the related relative position
        sensorD = 0.2347
        coorComcam = sensorD * np.array([-1, 0 ,1])

        # Generate the fields point x, y coordinates
        fieldX = np.kron(coorComcam, np.ones(nRow))
        fieldY = np.kron(np.ones(nCol), coorComcam)
        self.setFieldXYinDeg(fieldX, fieldY)

    def getZkFromOpd(self, opdFitsFile, znTerms=22, obscuration=0.61):
        """
        
        Get the wavefront error of optical path difference (OPD) in the basis of 
        annular Zernike polynomials.
        
        Arguments:
            opdFitsFile {[str]} -- OPD FITS file.
        
        Keyword Arguments:
            znTerms {int} -- Number of terms of annular Zk (z1-z22 by default). (default: {22})
            obscuration {float} -- Obscuration of annular Zernike polynomial. (default: {0.61})
        
        Returns:
            [ndarray] -- Annular Zernike polynomials. For PhoSim OPD, the unit is um.
            [ndarray] -- OPD map.
            [ndarray] -- Meshgrid x in OPD map.
            [ndarray] -- Meshgrid y in OPD map.
        
        Raises:
            RuntimeError -- The x, y dimensions of OPD are different.
        """

        # Get the OPD data (PhoSim OPD unit: um)
        opd = fits.getdata(opdFitsFile)

        # Check the x, y dimensions of OPD are the same
        if (np.unique(opd.shape).size != 1):
            raise RuntimeError("The x, y dimensions of OPD are different.")

        # x-, y-coordinate in the OPD image
        opdSize = opd.shape[0]
        opdGrid1d = np.linspace(-1, 1, opdSize)
        opdx, opdy = np.meshgrid(opdGrid1d, opdGrid1d)

        # Fit the OPD map with Zk and write into the file
        idx = (opd != 0)
        zk = ZernikeAnnularFit(opd[idx], opdx[idx], opdy[idx], znTerms, obscuration)

        return zk, opd, opdx, opdy

    def rmPTTfromOPD(self, opdFitsFile):
        """
        
        Remove the afftection of piston (z1), x-tilt (z2), and y-tilt (z3) from the optical 
        map difference (OPD) map.
        
        Arguments:
            opdFitsFile {[str]} -- OPD FITS file.
        
        Returns:
            [ndarray] -- OPD map after removing the affection of z1-z3.
            [ndarray] -- Meshgrid x in OPD map.
            [ndarray] -- Meshgrid y in OPD map.  
        """
        
        # Do the spherical Zernike fitting for the OPD map
        # Only fit the first three terms (z1-z3): piston, x-tilt, y-tilt
        zk, opd, opdx, opdy = self.getZkFromOpd(opdFitsFile, znTerms=3, obscuration=0)

        # Find the index that the value of OPD is not 0
        idx = (opd != 0)

        # Remove the PTT
        opd[idx] -= ZernikeEval(zk, opdx[idx], opdy[idx])

        return opd, opdx, opdy

    def addFieldXYbyCamPos(self, sensorName, xInpixel, yInPixel, folderPath2FocalPlane=None):
        """
        
        Add the new field X, Y in degree by the camera pixel positions.
        
        Arguments:
            sensorName {[str]} -- Canera sensor name (e.g. "R22_S11", "R40_S02_C0").
            xInpixel {[float]} -- Pixel x on camera coordinate.
            yInPixel {[float]} -- Pixel y on camera coordinate.
        
        Keyword Arguments:
            folderPath2FocalPlane {[str]} -- Path to the directory of focal plane data 
                                            ("focalplanelayout.txt") (default: {None})
        """

        # Get the focal plane data and set the sensor name
        sourProc = SourceProcessor()
        sourProc.config(sensorName=sensorName, folderPath2FocalPlane=folderPath2FocalPlane)

        # Do the coordinate transformation
        fieldXInDegree, fieldYInDegree = sourProc.camXYtoFieldXY(xInpixel, yInPixel)

        # Add to listed field x, y
        self.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

    def calcPSSN(self, opdFitsFile, wavelengthInUm, zen=0, debugLevel=0):
        """
        
        Calculate the normalized point source sensitivity (PSSN) based on the optical path 
        difference (OPD) map.
        
        Arguments:
            opdFitsFile {[str]} -- OPD FITS file.
            wavelengthInUm {[float]} -- Wavelength in microns.
        
        Keyword Arguments:
            zen {float} -- Telescope zenith angle in degree. (default: {0})
            debugLevel {int} -- Debug level. The higher value gives more information. (default: {0})
        
        Returns:
            [float] -- Calculated PSSN.
        """
        
        # Before calc_pssn,
        # (1) Remove PTT (piston, x-tilt, y-tilt),
        # (2) Make sure outside of pupil are all zeros
        opdRmPTT = self.rmPTTfromOPD(opdFitsFile)[0]

        # Calculate the normalized point source sensitivity (PSSN)
        pssn = calc_pssn(opdRmPTT, wavelengthInUm, zen=zen, debugLevel=debugLevel)

        return pssn

    def calcFWHMeff(self, pssn):
        """
        
        Calculate the effective full width at half maximum (FWHM).
        
        Arguments:
            pssn {[float]} -- Normalized point source sensitivity (PSSN).
        
        Returns:
            [float] -- Effective FWHM.
        """
        
        # FWHMeff_sys = FWHMeff_atm * sqrt(1/PSSN - 1). FWHMeff_atm = 0.6 arcsec. 
        # Another correction factor (eta = 1.086) is used to account for the difference between the 
        # simple RSS and the more proper convolution.
        # Follow page 7 (section 7.2 FWHM) in document-17242 for more information.
        eta = 1.086
        FWHMatm = 0.6
        FWHMeff = eta*FWHMatm*np.sqrt(1/pssn - 1)

        return FWHMeff

    def calcDm5(self, pssn):
        """
        
        Calculate the loss of limiting depth (dm5).
        
        Arguments:
            pssn {[float]} -- Normalized point source sensitivity (PSSN).
        
        Returns:
            [float] -- dm5.
        """
        
        # Calculate dm5 (the loss of limiting depth)
        # Check eq. (4.1) in page 4 in document-17242 for more information.
        dm5 = -1.25 * np.log10(pssn)

        return dm5

    def calcEllip(self, opdFitsFile, wavelengthInUm, zen=0, debugLevel=0):
        """
        
        Calculate the ellipticity.
        
        Arguments:
            opdFitsFile {[str]} -- OPD FITS file.
            wavelengthInUm {[float]} -- Wavelength in microns.
        
        Keyword Arguments:
            zen {float} -- Telescope zenith angle in degree. (default: {0})
            debugLevel {int} -- Debug level. The higher value gives more information. (default: {0})
        
        Returns:
            [float] -- Ellipticity.
        """

        # Remove the affection of piston (z1), x-tilt (z2), and y-tilt (z3) from OPD map.
        opdRmPTT = self.rmPTTfromOPD(opdFitsFile)[0]

        # Calculate the ellipticity
        elli = psf2eAtmW(opdRmPTT, wavelengthInUm, zen=zen, debugLevel=debugLevel)[0]

        return elli

    def calcGQvalue(self, valueList):
        """
        
        Calculate the value on Gaussian quadrature (GQ).
        
        Arguments:
            valueList {[list/ ndarray]} -- List of value (PSSN, effective FWHM, dm5, ellipticity).
        
        Returns:
            [float] -- GQ value.
        
        Raises:
            RuntimeError -- Length of weighting ratio != length of value list.
        """

        # Check the lengths of weighting ratio and value list are the same
        if (len(self.wt) != len(valueList)):
            raise RuntimeError("Length of weighting ratio != length of value list.")

        # Calculate the effective value on Gaussain quardure plane
        GQvalue = np.sum(self.wt*valueList)

        return GQvalue

if __name__ == "__main__":

    metr = OpdMetrology()

    # Calculate Zk based on OPD
    opdFitsFile = "/Users/Wolf/Documents/aosOutput/image/sim6/iter0/sim6_iter0_opd0.fits.gz"
    zk = metr.getZkFromOpd(opdFitsFile)[0]

    # Folder path of focal plane data
    folderPath2FocalPlane = "/Users/Wolf/Documents/bitbucket/phosim_syseng2/data/lsst"
    sensorName = "R13_S11"
    xInpixel = 2000
    yInPixel = 2036
    metr.addFieldXYbyCamPos(sensorName, xInpixel, yInPixel, folderPath2FocalPlane=folderPath2FocalPlane)

    # Remove the PTT
    opdRmPTT, opdx, opdy = metr.rmPTTfromOPD(opdFitsFile)

    # Calculate PSSN
    wavelengthInUm = 0.5
    pssn = metr.calcPSSN(opdFitsFile, wavelengthInUm)
    print(pssn)

    # Calculate the GQ PSSN
    metr.setDefaultLsstGQ()
    pssnList = []
    for ii in range(31):
        opdFitsFile = "/Users/Wolf/Documents/aosOutput/image/sim6/iter0/sim6_iter0_opd%d.fits.gz" % ii
        pssn = metr.calcPSSN(opdFitsFile, wavelengthInUm)
        pssnList.append(pssn)
    GQpssn = metr.calcGQvalue(pssnList)
    print(GQpssn)

    # Calculate the effective FWHM
    FWHMeff = metr.calcFWHMeff(pssn)
    print(FWHMeff)

    # Calculate the dm5
    dm5 = metr.calcDm5(pssn)
    print(dm5)

    # Calculate the ellipticity
    elli = metr.calcEllip(opdFitsFile, wavelengthInUm)
    print(elli)
