import os, unittest
import numpy as np
from astropy.io import fits
import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lsst.ts.wep.cwfs.Tool import ZernikeAnnularFit, ZernikeEval
from lsst.ts.wep.SourceProcessor import SourceProcessor
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

    def getZkFromOpd(self, opdFitsFile=None, opdMap=None, znTerms=22, obscuration=0.61):
        """
        
        Get the wavefront error of optical path difference (OPD) in the basis of 
        annular Zernike polynomials.
                
        Keyword Arguments:
            opdFitsFile {[str]} -- OPD FITS file. (default: {None})
            opdMap {[ndarray]} -- OPD map data. (default: {None})
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
        if (opdFitsFile is not None):
            opd = fits.getdata(opdFitsFile)
        elif (opdMap is not None):
            opd = opdMap.copy()

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

    def rmPTTfromOPD(self, opdFitsFile=None, opdMap=None):
        """
        
        Remove the afftection of piston (z1), x-tilt (z2), and y-tilt (z3) from the optical 
        map difference (OPD) map.
        
        Keyword Arguments:
            opdFitsFile {[str]} -- OPD FITS file. (default: {None})
            opdMap {[ndarray]} -- OPD map data. (default: {None})
        
        Returns:
            [ndarray] -- OPD map after removing the affection of z1-z3.
            [ndarray] -- Meshgrid x in OPD map.
            [ndarray] -- Meshgrid y in OPD map.  
        """
        
        # Do the spherical Zernike fitting for the OPD map
        # Only fit the first three terms (z1-z3): piston, x-tilt, y-tilt
        zk, opd, opdx, opdy = self.getZkFromOpd(opdFitsFile=opdFitsFile, opdMap=opdMap, 
                                                znTerms=3, obscuration=0)

        # Find the index that the value of OPD is not 0
        idx = (opd != 0)

        # Remove the PTT
        opd[idx] -= ZernikeEval(zk, opdx[idx], opdy[idx])

        return opd, opdx, opdy

    def addFieldXYbyCamPos(self, sensorName, xInpixel, yInPixel, folderPath2FocalPlane):
        """
        
        Add the new field X, Y in degree by the camera pixel positions.
        
        Arguments:
            sensorName {[str]} -- Canera sensor name (e.g. "R22_S11", "R40_S02_C0").
            xInpixel {[float]} -- Pixel x on camera coordinate.
            yInPixel {[float]} -- Pixel y on camera coordinate.
            folderPath2FocalPlane {[str]} -- Path to the directory of focal plane data 
                                            ("focalplanelayout.txt").
        """

        # Get the focal plane data and set the sensor name
        sourProc = SourceProcessor()
        sourProc.config(sensorName=sensorName, folderPath2FocalPlane=folderPath2FocalPlane)

        # Do the coordinate transformation
        fieldXInDegree, fieldYInDegree = sourProc.camXYtoFieldXY(xInpixel, yInPixel)

        # Add to listed field x, y
        self.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

    def calcPSSN(self, wavelengthInUm, opdFitsFile=None, opdMap=None, zen=0, debugLevel=0):
        """
        
        Calculate the normalized point source sensitivity (PSSN) based on the optical path 
        difference (OPD) map.
        
        Arguments:
            wavelengthInUm {[float]} -- Wavelength in microns.
        
        Keyword Arguments:
            opdFitsFile {[str]} -- OPD FITS file. (default: {None})
            opdMap {[ndarray]} -- OPD map data. (default: {None})
            zen {float} -- Telescope zenith angle in degree. (default: {0})
            debugLevel {int} -- Debug level. The higher value gives more information. (default: {0})
        
        Returns:
            [float] -- Calculated PSSN.
        """
        
        # Before calc_pssn,
        # (1) Remove PTT (piston, x-tilt, y-tilt),
        # (2) Make sure outside of pupil are all zeros
        opdRmPTT = self.rmPTTfromOPD(opdFitsFile=opdFitsFile, opdMap=opdMap)[0]

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

    def calcEllip(self, wavelengthInUm, opdFitsFile=None, opdMap=None, zen=0, debugLevel=0):
        """
        
        Calculate the ellipticity.
        
        Arguments:
            wavelengthInUm {[float]} -- Wavelength in microns.
        
        Keyword Arguments:
            opdFitsFile {[str]} -- OPD FITS file. (default: {None})
            opdMap {[ndarray]} -- OPD map data. (default: {None})
            zen {float} -- Telescope zenith angle in degree. (default: {0})
            debugLevel {int} -- Debug level. The higher value gives more information. (default: {0})
        
        Returns:
            [float] -- Ellipticity.
        """

        # Remove the affection of piston (z1), x-tilt (z2), and y-tilt (z3) from OPD map.
        opdRmPTT = self.rmPTTfromOPD(opdFitsFile=opdFitsFile, opdMap=opdMap)[0]

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

    def showFieldMap(self, folderPath2FocalPlane=None, saveToFilePath=None, dpi=None, pixel2Arcsec=0.2):
        """
        
        Show the field map in degree.
        
        Keyword Arguments:
            folderPath2FocalPlane {[str]} -- Folder directory to focal plane file. 
                                            (default: {None})
            saveToFilePath {str} -- File path to save the figure. (default: {None})
            dpi {int} -- The resolution in dots per inch. (default: {None})
            pixel2Arcsec {float} -- Pixel to arcsec. (default: {0.2})
        """

        # Declare the figure
        plt.figure()
        
        # Get the focal plane information
        if (folderPath2FocalPlane is not None):
            sourProc = SourceProcessor()
            sourProc.config(folderPath2FocalPlane=folderPath2FocalPlane)

            # Plot the CCD boundary
            for sensorName in sourProc.sensorDimList.keys():

                # Get the CCD corner field points in degree
                pointXinDeg, pointYinDeg = self.__getCCDBoundInDeg(sourProc, sensorName, 
                                                                pixel2Arcsec=pixel2Arcsec)

                # Plot the boundary
                pointXinDeg.append(pointXinDeg[0])
                pointYinDeg.append(pointYinDeg[0])
                plt.plot(pointXinDeg, pointYinDeg, "b")

        # Plot the field X, Y position
        plt.plot(self.fieldX, self.fieldY, "ro")

        # Plot the 3.5 degree circle
        circle = plt.Circle((0, 0), 1.75, color="g", fill=False)
        plt.gcf().gca().add_artist(circle)

        # Do the labeling
        plt.xlabel("Field X (deg)")
        plt.ylabel("Field Y (deg)")

        # Label four corner rafts for the identification
        plt.text(-1.5, -1.5, "R00")
        plt.text(1.5, -1.5, "R40")
        plt.text(1.5, 1.5, "R44")
        plt.text(-1.5, 1.5, "R04")

        # Set the axis limit
        plt.xlim(-2, 2)
        plt.ylim(-2, 2)

        # Set the same scale
        plt.axis("equal")

        # Save the figure or not
        if (saveToFilePath is not None):
            plt.savefig(saveToFilePath, dpi=dpi)
            plt.close()
        else:
            plt.show()

    def __getCCDBoundInDeg(self, sourProc, sensorName, pixel2Arcsec=0.2):
        """
        
        Get the CCD four corners in degree in counter clockwise direction.
        
        Arguments:
            sourProc {[SourceProcessor]} -- SourceProcessor object.
            sensorName {[str]} -- Sensor name
        
        Keyword Arguments:
            pixel2Arcsec {float} -- Pixel to arcsec. (default: {0.2})
        
        Returns:
            [list] -- Corner points x in degree.
            [list] -- Corner points y in degree.
        """

        # Set the sensor on sourProc
        sourProc.config(sensorName=sensorName, pixel2Arcsec=pixel2Arcsec)

        # Get the dimension of CCD
        pixDimX, pixDimY = sourProc.sensorDimList[sensorName]

        # Define four corner points in the counter-clockwise direction
        pixelXlist = [0, pixDimX, pixDimX, 0]
        pixelYlist = [0, 0, pixDimY, pixDimY]

        # Transform from pixel to field degree
        fieldXinDegList = []
        fieldYinDeglist = []

        for ii in range(len(pixelXlist)):

            fieldX, fieldY = sourProc.camXYtoFieldXY(pixelXlist[ii], pixelYlist[ii])
            fieldXinDegList.append(fieldX)
            fieldYinDeglist.append(fieldY)

        return fieldXinDegList, fieldYinDeglist

class OpdMetrologyTest(unittest.TestCase):
    
    """
    Test functions in OpdMetrology.
    """

    def setUp(self):

        self.testDataDir = os.path.join("..", "testData", "testOpdFunc")

    def testFunc(self):

        # Instantiate the OpdMetrology object
        metr = OpdMetrology()

        wt = np.array([1])
        metr.setWeightingRatio(wt)
        self.assertEqual(metr.wt, wt)

        fieldXInDegree = 0.1
        fieldYInDegree = 0.2
        metr.setFieldXYinDeg(fieldXInDegree, fieldYInDegree)
        self.assertEqual(metr.fieldX, fieldXInDegree)
        self.assertEqual(metr.fieldY, fieldYInDegree)

        metr.addWeightingRatio(1)
        self.assertEqual(len(metr.wt), 2)

        metr.normalizeWeightingRatio()
        self.assertEqual(np.sum(metr.wt), 1)

        fieldXInDegree = 0.2
        fieldYInDegree = 0.3
        metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)
        self.assertEqual(len(metr.fieldX), 2)

        metr.setDefaultLsstGQ()
        self.assertEqual(len(metr.fieldX), 31)

        fieldWFSx, fieldWFSy = metr.getDefaultLsstWfsGQ()
        self.assertEqual(len(fieldWFSx), 4)

        metr.setDefaultComcamGQ()
        self.assertEqual(len(metr.fieldX), 9)

        opdFileName = "sim6_iter0_opd0.fits.gz"
        opdFilePath = os.path.join(self.testDataDir, opdFileName)
        zk = metr.getZkFromOpd(opdFitsFile=opdFilePath)[0]
        
        ansOpdFileName = "sim6_iter0_opd.zer"
        ansOpdFilePath = os.path.join(self.testDataDir, ansOpdFileName)
        allOpdAns = np.loadtxt(ansOpdFilePath)
        self.assertLess(np.sum(np.abs(zk-allOpdAns[0,:])), 1e-10)

        opdRmPTT, opdx, opdy = metr.rmPTTfromOPD(opdFitsFile=opdFilePath)
        zkRmPTT = metr.getZkFromOpd(opdMap=opdRmPTT)[0]
        self.assertLess(np.sum(np.abs(zkRmPTT[0:3])), 5e-2)

        sensorName = "R22_S11"
        xInpixel = 4000
        yInPixel = 4072
        metr.addFieldXYbyCamPos(sensorName, xInpixel, yInPixel, self.testDataDir)

        ansFieldXinDeg = 2000*0.2/3600
        ansFieldYinDeg = 2036*0.2/3600
        self.assertAlmostEqual((metr.fieldX[-1], metr.fieldY[-1]), (ansFieldXinDeg, ansFieldYinDeg))

        wavelengthInUm = 0.5
        pssn = metr.calcPSSN(wavelengthInUm, opdFitsFile=opdFilePath)

        ansAllDataFileName = "sim6_iter0_PSSN.txt"
        ansAllDataFilePath = os.path.join(self.testDataDir, ansAllDataFileName)
        allData = np.loadtxt(ansAllDataFilePath)
        self.assertAlmostEqual(pssn, allData[0,0])

        fwhm = metr.calcFWHMeff(pssn)
        self.assertAlmostEqual(fwhm, allData[1,0])

        dm5 = metr.calcDm5(pssn)
        self.assertAlmostEqual(dm5, allData[2,0])

        elli = metr.calcEllip(wavelengthInUm, opdFitsFile=opdFilePath)

        ansElliFileName = "sim6_iter0_elli.txt"
        ansElliFilePath = os.path.join(self.testDataDir, ansElliFileName)
        allElli = np.loadtxt(ansElliFilePath)
        self.assertAlmostEqual(elli, allElli[0])

        metr.setDefaultLsstGQ()
        valueList = allData[0, 0:31]
        GQvalue = metr.calcGQvalue(valueList)
        self.assertAlmostEqual(GQvalue, allData[0,-1])

        saveToFilePath = os.path.join("..", "outputImg", "fieldMap.png")
        metr.showFieldMap(folderPath2FocalPlane=self.testDataDir, saveToFilePath=saveToFilePath)
        os.remove(saveToFilePath)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
