import os
import numpy as np
from scipy.interpolate import Rbf

from lsst.ts.phosim.MirrorSim import MirrorSim
from lsst.ts.phosim.Utility import phosim2ZemaxCoorTrans
from lsst.ts.phosim.PlotUtil import plotResMap
from lsst.ts.wep.cwfs.Tool import ZernikeAnnularFit, ZernikeAnnularEval


class M1M3Sim(MirrorSim):
    
    def __init__(self, mirrorDataDir=""):
        """Initiate the M1M3 simulator class.

        Parameters
        ----------
        mirrorDataDir : {str}, optional
            Mirror data directory. (the default is "".)
        """

        # Inner and outer radius of M1 mirror in m
        R1i = 2.558
        R1 = 4.180

        # Inner and outer radius of M3 mirror in m
        R3i = 0.550
        R3 = 2.508

        super(M1M3Sim, self).__init__((R1i, R3i), (R1, R3),
                                      mirrorDataDir=mirrorDataDir)

    def getActForce(self, actForceFileName="M1M3_1um_156_force.DAT"):
        """Get the mirror actuator forces in N.

        Parameters
        ----------
        actForceFileName : str, optional
            Actuator force file name. (the default is "M1M3_1um_156_force.DAT".)

        Returns
        -------
        numpy.ndarray
            Actuator forces in N.
        """

        return super(M1M3Sim, self).getActForce(actForceFileName)

    def getPrintthz(self, zAngleInRadian, preCompElevInRadian=0,
                    FEAfileName="", FEAzenFileName="M1M3_dxdydz_zenith.txt",
                    FEAhorFileName="M1M3_dxdydz_horizon.txt",
                    gridFileName="M1M3_1um_156_grid.DAT"):
        """Get the mirror print in m along z direction in specific zenith
        angle.

        FEA: Finite element analysis.

        Parameters
        ----------
        zAngleInRadian : float
            Zenith angle in radian.
        preCompElevInRadian : float, optional
            Pre-compensation elevation angle in radian. (the default is 0.)
        FEAfileName : str, optional
            FEA model data file name. (the default is "".)
        FEAzenFileName : str, optional
            FEA model data file name in zenith angle. (the default is
            "M1M3_dxdydz_zenith.txt".)
        FEAhorFileName : str, optional
            FEA model data file name in horizontal angle. (the default is
            "M1M3_dxdydz_horizon.txt".)
        gridFileName : str, optional
            File name of bending mode data. (the default is
            "M1M3_1um_156_grid.DAT".)
        
        Returns
        -------
        numpy.ndarray
            Corrected projection in m along z direction.
        """

        # Data needed to determine gravitational print through
        data = self.getMirrorData(FEAzenFileName)
        zdx = data[:, 0]
        zdy = data[:, 1]
        zdz = data[:, 2]

        data = self.getMirrorData(FEAhorFileName)
        hdx = data[:, 0]
        hdy = data[:, 1]
        hdz = data[:, 2]

        # Do the M1M3 gravitational correction.
        # Map the changes of dx, dy, and dz on a plane for certain zenith angle
        printthxInM = zdx * np.cos(zAngleInRadian) + \
                      hdx * np.sin(zAngleInRadian)
        printthyInM = zdy * np.cos(zAngleInRadian) + \
                      hdy * np.sin(zAngleInRadian)
        printthzInM = zdz * np.cos(zAngleInRadian) + \
                      hdz * np.sin(zAngleInRadian)

        # Get the bending mode information
        idx1, idx3, bx, by, bz = self._getMirCoor(gridFileName)

        # Calcualte the mirror ideal shape
        zRef = self._calcIdealShape(bx*1000, by*1000, idx1, idx3)/1000

        # Calcualte the mirror ideal shape with the displacement
        zpRef = self._calcIdealShape((bx + printthxInM)*1000,
                                     (by + printthyInM)*1000,
                                     idx1, idx3)/1000

        # Convert printthz into surface sag to get the estimated wavefront
        # error.
        # Do the zenith angle correction by the linear approximation with the
        # ideal shape.
        printthzInM = printthzInM - (zpRef - zRef)

        # Normalize the coordinate
        Ri = self.getInnerRinM()[0]
        R = self.getOuterRinM()[0]

        normX = bx/R
        normY = by/R
        obs = Ri/R

        # Fit the annular Zernike polynomials z0-z2 (piton, x-tilt, y-tilt)
        zc = ZernikeAnnularFit(printthzInM, normX, normY, 3, obs)

        # Do the estimated wavefront error correction for the mirror projection
        printthzInM -= ZernikeAnnularEval(zc, normX, normY, obs)

        return printthzInM

    def _getMirCoor(self, gridFileName):
        """Get the mirror coordinate and node.

        Parameters
        ----------
        gridFileName : str
            File name of bending mode data.

        Returns
        -------
        numpy.ndarray[int]
            M1 node.
        numpy.ndarray[int]
            M3 node.
        numpy.ndarray
            x coordinate.
        numpy.ndarray
            y coordinate.
        numpy.ndarray
            z coordinate.
        """

        # Get the bending mode information
        data = self.getMirrorData(gridFileName)
        
        nodeID = data[:, 0].astype("int")
        nodeM1 = (nodeID == 1)
        nodeM3 = (nodeID == 3)
        
        bx = data[:, 1]
        by = data[:, 2]
        bz = data[:, 3:]

        return nodeM1, nodeM3, bx, by, bz

    def _calcIdealShape(self, xInMm, yInMm, idxM1, idxM3, dr1=0, dr3=0, dk1=0,
                        dk3=0):
        """Calculate the ideal shape of mirror along z direction.
        
        This is described by a series of cylindrically-symmetric aspheric
        surfaces.
        
        Parameters
        ----------
        xInMm : numpy.ndarray
            Coordinate x in 1D array in mm.
        yInMm : numpy.ndarray
            Coordinate y in 1D array in mm.
        idxM1 : numpy.ndarray [int]
            M1 node.
        idxM3 : numpy.ndarray [int]
            M3 node.
        dr1 : float, optional
            Displacement of r in mirror 1. (the default is 0.)
        dr3 : float, optional
            Displacement of r in mirror 3. (the default is 0.)
        dk1 : float, optional
            Displacement of kappa (k) in mirror 1. (the default is 0.)
        dk3 : float, optional
            Displacement of kappa (k) in mirror 3. (the default is 0.)

        Returns
        -------
        numpy.ndarray
            Ideal mirror surface along z direction.
        """

        # M1 optical design
        r1 = -1.9835e4
        k1 = -1.215
        alpha1 = np.zeros((8, 1))
        alpha1[2] = 1.38e-24

        # M3 optical design
        r3 = -8344.5
        k3 = 0.155
        alpha3 = np.zeros((8, 1))
        alpha3[2] = -4.5e-22
        alpha3[3] = -8.15e-30

        # Get the dimension of input xInMm, yInMm
        nr = xInMm.shape
        mr = yInMm.shape
        if (nr != mr):
            print("In the ideal shape calculation, x is [%d] while y is [%d]."
                  % (nr, mr))
            sys.exit()

        # Calculation the curvature (c) and conic constant (kappa)

        # Mirror 1 (M1)
        c1 = 1/(r1 + dr1)
        k1 = k1 + dk1

        # Mirror 3 (M3)
        c3 = 1/(r3 + dr3)
        k3 = k3 + dk3

        # Construct the curvature, kappa, and alpha matrixes for the ideal
        # shape calculation
        cMat = np.zeros(nr)
        cMat[idxM1] = c1
        cMat[idxM3] = c3

        kMat = np.zeros(nr)
        kMat[idxM1] = k1
        kMat[idxM3] = k3

        alphaMat = np.tile(np.zeros(nr), (8, 1))
        for ii in range(8):
            alphaMat[ii, idxM1] = alpha1[ii]
            alphaMat[ii, idxM3] = alpha3[ii]

        # Calculate the radius
        r2 = xInMm**2 + yInMm**2
        r = np.sqrt(r2)

        # Calculate the ideal surface

        # The optical elements of telescopes can often be described by a
        # series of cylindrically-symmetric aspheric surfaces:
        # z(r) = c * r^2/[ 1 + sqrt( 1-(1+k) * c^2 * r^2 ) ] + sum(ai * r^(2*i)) + sum(Aj * Zj)
        # where i = 1-8, j = 1-N

        z0 = cMat * r2 / (1 + np.sqrt(1 - (1 + kMat) * cMat**2 * r2))
        for ii in range(8):
            z0 += alphaMat[ii, :] * r2**(ii+1)

        # M3 vertex offset from M1 vertex, values from Zemax model
        # M3voffset = (233.8 - 233.8 - 900 - 3910.701 - 1345.500 + 1725.701
        # + 3530.500 + 900 + 233.800)
        M3voffset = 233.8

        # Add the M3 offset (sum(Aj * Zj), j = 1 - N)
        z0[idxM3] = z0[idxM3] + M3voffset

        # In Zemax, z axis points from M1M3 to M2. the reversed direction
        # (z0>0) is needed. That means the direction of M2 to M1M3.
        return -z0

    def getTempCorr(self, M1M3TBulk, M1M3TxGrad, M1M3TyGrad, M1M3TzGrad,
                    M1M3TrGrad, FEAfileName="M1M3_thermal_FEA.txt",
                    gridFileName="M1M3_1um_156_grid.DAT"):
        """Get the mirror print correction along z direction for certain
        temperature gradient.

        FEA: Finite element analysis.

        Parameters
        ----------
        M1M3TBulk : float
            Bulk temperature in degree C (+/-2sigma spans +/-0.8C).
        M1M3TxGrad : float
            Temperature gradient along x direction in degree C (+/-2sigma
            spans 0.4C).
        M1M3TyGrad : float
            Temperature gradient along y direction in degree C (+/-2sigma
            spans 0.4C).
        M1M3TzGrad : float
            Temperature gradient along z direction in degree C (+/-2sigma
            spans 0.1C).
        M1M3TrGrad : float
            Temperature gradient along r direction in degree C (+/-2sigma
            spans 0.1C).
        FEAfileName : str, optional
            FEA model data file name.  (the default is "M1M3_thermal_FEA.txt".)
        gridFileName : str, optional
            File name of bending mode data. (the default is
            "M1M3_1um_156_grid.DAT".)

        Returns
        -------
        numpy.ndarray
            Corrected projection in um along z direction.
        """

        # Data needed to determine thermal deformation
        data = self.getMirrorData(FEAfileName, skiprows=1)

        # These are the normalized coordinates

        # In the original XLS file (thermal FEA model), max(x)=164.6060 in,
        # while 4.18m = 164.5669 in for real mirror. These two numbers do
        # not match.
        # The data here has normalized the dimension (x, y) already.
        # n.b. these may not have been normalized correctly, b/c max(tx)=1.0
        tx = data[:, 0]
        ty = data[:, 1]

        # Below are in M1M3 coordinate system, and in micron

        # Do the fitting in the normalized coordinate
        bx, by = self._getMirCoor(gridFileName)[2:4]
        R = self.RinM[0]
        normX = bx/R
        normY = by/R

        # Fit the bulk
        tbdz = self._fitData(tx, ty, data[:, 2], normX, normY)

        # Fit the x-grad
        txdz = self._fitData(tx, ty, data[:, 3], normX, normY)

        # Fit the y-grad
        tydz = self._fitData(tx, ty, data[:, 4], normX, normY)

        # Fit the z-grad
        tzdz = self._fitData(tx, ty, data[:, 5], normX, normY)

        # Fit the r-gradß
        trdz = self._fitData(tx, ty, data[:, 6], normX, normY)

        # Get the temprature correction
        tempCorrInUm = M1M3TBulk*tbdz + M1M3TxGrad*txdz + M1M3TyGrad*tydz + \
                       M1M3TzGrad*tzdz + M1M3TrGrad*trdz

        return tempCorrInUm

    def _fitData(self, dataX, dataY, data, x, y):
        """Fit the data by radial basis function.

        Parameters
        ----------
        dataX : numpy.ndarray
            Data x.
        dataY : numpy.ndarray
            Data y.
        data : numpy.ndarray
            Data to fit.
        x : numpy.ndarray
            x coordinate.
        y : numpy.ndarray
            y coordinate.

        Returns
        -------
        numpy.ndarray
            Fitted data.
        """

        # Construct the fitting model
        rbfi = Rbf(dataX, dataY, data)

        # Return the fitted data
        return rbfi(x, y)

    def getMirrorResInMmInZemax(self, gridFileName="M1M3_1um_156_grid.DAT",
                                numTerms=28, writeZcInMnToFilePath=None):
        """Get the residue of surface (mirror print along z-axis) in mm under
        the Zemax coordinate.

        This value is after the fitting with spherical Zernike polynomials
        (zk).

        Parameters
        ----------
        gridFileName : str, optional
            File name of bending mode data. (the default is
            "M1M3_1um_156_grid.DAT".)
        numTerms : {number}, optional
            Number of Zernike terms to fit. (the default is 28.)
        writeZcInMnToFilePath : str, optional
            File path to write the fitted zk in mm. (the default is None.)

        Returns
        ------
        numpy.ndarray
            Fitted residue in mm after removing the fitted zk terms in Zemax
            coordinate.
        numpy.ndarray
            X position in mm in Zemax coordinate.
        numpy.ndarray
            Y position in mm in Zemax coordinate.
        numpy.ndarray
            Fitted zk in mm in Zemax coordinate.
        """

        # Get the bending mode information
        bx, by, bz = self._getMirCoor(gridFileName)[2:5]

        # Transform the M1M3 coordinate to Zemax coordinate
        bxInZemax, byInZemax, surfInZemax = phosim2ZemaxCoorTrans(
                                                bx, by, self.getSurfAlongZ())

        # Get the mirror residue and zk in um
        RinM = self.getOuterRinM()[0]
        resInUmInZemax, zcInUmInZemax = self._getMirrorResInNormalizedCoor(
                    surfInZemax, bxInZemax/RinM, byInZemax/RinM, numTerms)

        # Change the unit to mm
        resInMmInZemax = resInUmInZemax * 1e-3
        bxInMmInZemax = bxInZemax * 1e3
        byInMmInZemax = byInZemax * 1e3
        zcInMmInZemax = zcInUmInZemax * 1e-3

        # Save the file of fitted Zk
        if (writeZcInMnToFilePath is not None):
            np.savetxt(writeZcInMnToFilePath, zcInMmInZemax)

        return resInMmInZemax, bxInMmInZemax, byInMmInZemax, zcInMmInZemax

    def writeMirZkAndGridResInZemax(self, resFile=[], surfaceGridN=200,
                                    gridFileName="M1M3_1um_156_grid.DAT",
                                    numTerms=28, writeZcInMnToFilePath=None):
        """
        
        Write the grid residue in mm of mirror surface after the fitting with Zk under the Zemax
        coordinate.
        
        Keyword Arguments:
            resFile {[list]} -- File path to save the grid surface residue map ([M1filePath, M3filePath]). 
                                (default: {[]]})
            surfaceGridN {int} -- Surface grid number. (default: {200})
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
            numTerms {int} -- Number of Zernike terms to fit. (default: {28})
            writeZcInMnToFilePath {[str]} -- File path to write the fitted zk in mm. (default: {None})
        
        Returns:
            [str] -- Grid residue map related data of M1.
            [str] -- Grid residue map related data of M3.
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = \
            self.getMirrorResInMmInZemax(
                        gridFileName=gridFileName, numTerms=numTerms,
                        writeZcInMnToFilePath=writeZcInMnToFilePath)[0:3]

        # Get the mirror node
        idx1, idx3 = self._getMirCoor(gridFileName)[0:2]

        # Grid sample map for M1 and M3
        for ii, idx in zip((0, 1), (idx1, idx3)):

            # Change the unit from m to mm
            innerRinMm = self.getInnerRinM()[ii] * 1e3
            outerRinMm = self.getOuterRinM()[ii] * 1e3

            # Get the residue map used in Zemax
            # Content header: (NUM_X_PIXELS, NUM_Y_PIXELS, delta x, delta y)
            # Content: (z, dx, dy, dxdy)
            content = self._gridSampInMnInZemax(
                                resInMmInZemax[idx], bxInMmInZemax[idx], 
                                byInMmInZemax[idx], innerRinMm, outerRinMm, 
                                surfaceGridN, surfaceGridN,
                                resFile=resFile[ii])
            if (ii == 0):
                contentM1 = content
            elif (ii == 1):
                contentM3 = content

        return contentM1, contentM3

    def showMirResMap(self, gridFileName="M1M3_1um_156_grid.DAT", numTerms=28,
                      resFile=[], writeToResMapFilePath=[]):
        """
        
        Show the mirror residue map.
        
        Keyword Arguments:
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
            numTerms {int} -- Number of Zernike terms to fit. (default: {28})
            resFile {list} -- File path of the grid surface residue map. (default: {[]})
            writeToResMapFilePath {list} -- File path to save the residue map. (default: {[]})
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = \
            self.getMirrorResInMmInZemax(gridFileName=gridFileName,
                                         numTerms=numTerms)[0:3]

        # Get the mirror node
        idx1, idx3 = self._getMirCoor(gridFileName)[0:2]

        # Show the mirror maps
        RinMtuple = self.getOuterRinM()

        for ii, RinM, idx in zip((0, 1), RinMtuple, (idx1, idx3)):
            outerRinMm = RinM * 1e3
            plotResMap(resInMmInZemax[idx], bxInMmInZemax[idx],
                       byInMmInZemax[idx], outerRinMm, resFile=resFile[ii],
                       writeToResMapFilePath=writeToResMapFilePath[ii])

    def genMirSurfRandErr(self, zAngleInRadian, LUTfileName="M1M3_LUT.txt", 
                          forceZenFileName="M1M3_force_zenith.txt", 
                          forceHorFileName="M1M3_force_horizon.txt", 
                          forceInflFileName="M1M3_influence_256.txt", 
                          M1M3ForceError=0.05, nzActuator=156, seedNum=0):
        """Generate the mirror surface random error.

        LUT: Loop-up table.

        Parameters
        ----------
        zAngleInRadian : float
            Zenith angle in radian.
        LUTfileName : str, optional
            LUT file name. (the default is "M1M3_LUT.txt".)
        forceZenFileName : str, optional
            File name of actuator forces along zenith direction. (the default
            is "M1M3_force_zenith.txt".)
        forceHorFileName : str, optional
            File name of actuator forces along horizon direction. (the default
            is "M1M3_force_horizon.txt".)
        forceInflFileName : str, optional
            Influence matrix of actuator forces. (the default is
            "M1M3_influence_256.txt".)
        M1M3ForceError : float, optional
            Ratio of actuator force error. (the default is 0.05.)
        nzActuator : int, optional
            Number of actuator along z direction. (the default is 156.)
        seedNum : int, optional
            Random seed number. (the default is 0.)
        
        Returns
        -------
        numpy.ndarray
            Generated mirror surface random error in m.
        """

        # Get the actuator forces in N of M1M3 based on the look-up table (LUT)
        zangleInDeg = np.rad2deg(zAngleInRadian)
        LUTforce = self.getLUTforce(zangleInDeg, LUTfileName)

        # Add 5% force error (self.M1M3ForceError). This is for iteration 0
        # only.
        # This means from -5% to +5% of original actuator's force.
        np.random.seed(int(seedNum))
        nActuator = len(LUTforce)
        myu = (1 + 2*(np.random.rand(nActuator) - 0.5)*M1M3ForceError)*LUTforce

        # Balance forces along z-axis
        # This statement is intentionally to make the force balance.
        nzActuator = int(nzActuator)
        myu[nzActuator-1] = np.sum(LUTforce[:nzActuator]) - \
                            np.sum(myu[:nzActuator-1])

        # Balance forces along y-axis
        # This statement is intentionally to make the force balance.
        myu[nActuator-1] = np.sum(LUTforce[nzActuator:]) - \
                           np.sum(myu[nzActuator:-1])

        # Get the net force along the z-axis
        zf = self.getMirrorData(forceZenFileName)
        hf = self.getMirrorData(forceHorFileName)
        u0 = zf*np.cos(zAngleInRadian) + hf*np.sin(zAngleInRadian)

        # Calculate the random surface
        G = self.getMirrorData(forceInflFileName)
        randSurfInM = G.dot(myu - u0)

        return randSurfInM


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
    