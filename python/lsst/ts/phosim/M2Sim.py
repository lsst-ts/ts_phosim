import os
import numpy as np

from lsst.ts.phosim.MirrorSim import MirrorSim
from lsst.ts.phosim.Utility import phosim2ZemaxCoorTrans
from lsst.ts.phosim.PlotUtil import plotResMap


class M2Sim(MirrorSim):

    def __init__(self, mirrorDataDir=""):
        """Initiate the M2 simulator class.

        Parameters
        ----------
        mirrorDataDir : str, optional
            Mirror data directory. (the default is "".)
        """

        # Inner and outer radius of M2 mirror in m
        Ri = 0.9
        R = 1.710

        super(M2Sim, self).__init__(Ri, R, mirrorDataDir=mirrorDataDir)

    def getActForce(self, actForceFileName="M2_1um_force.DAT"):
        """Get the mirror actuator forces in N.

        Parameters
        ----------
        actForceFileName : str, optional
            Actuator force file name. (the default is "M2_1um_force.DAT".)
        
        Returns
        -------
        numpy.ndarray
            Actuator forces in N.
        """

        return super(M2Sim, self).getActForce(actForceFileName)

    def getPrintthz(self, zAngleInRadian, preCompElevInRadian=0,
                    FEAfileName="M2_GT_FEA.txt", FEAzenFileName="",
                    FEAhorFileName="", gridFileName=""):
        """Get the mirror print in um along z direction in specific zenith
        angle.

        FEA: Finite element analysis.

        Parameters
        ----------
        zAngleInRadian : float
            Zenith angle in radian.
        preCompElevInRadian : float, optional
            Pre-compensation elevation angle in radian. (the default is 0.)
        FEAfileName : str, optional
            FEA model data file name. (the default is "M2_GT_FEA.txt".)
        FEAzenFileName : str, optional
            FEA model data file name in zenith angle. (the default is "".)
        FEAhorFileName : str, optional
            FEA model data file name in horizontal angle. (the default is "".)
        gridFileName : str, optional
            File name of bending mode data. (the default is "".)

        Returns
        ------
        numpy.ndarray
            Corrected projection in um along z direction.
        """

        # Read the FEA file
        data = self.getMirrorData(FEAfileName, skiprows=1)

        # Zenith direction in um
        zdz = data[:, 2]

        # Horizon direction in um
        hdz = data[:, 3]

        # Do the M2 gravitational correction.
        # Map the changes of dz on a plane for certain zenith angle
        printthzInUm = zdz * np.cos(zAngleInRadian) + \
                       hdz * np.sin(zAngleInRadian)

        # Do the pre-compensation elevation angle correction
        printthzInUm -= zdz * np.cos(preCompElevInRadian) + \
                        hdz * np.sin(preCompElevInRadian)

        return printthzInUm

    def getTempCorr(self, M2TzGrad, M2TrGrad, FEAfileName="M2_GT_FEA.txt"):
        """Get the mirror print correction along z direction for certain
        temperature gradient.

        FEA: Finite element analysis.

        Parameters
        ----------
        M2TzGrad : float
            Temperature gradient along z direction in degree C. (+/-2sigma
            spans 1C).
        M2TrGrad : float
            Temperature gradient along r direction in degree C. (+/-2sigma
            spans 1C).
        FEAfileName : str, optional
            FEA model data file name. (the default is "M2_GT_FEA.txt".)
        
        Returns
        -------
        numpy.ndarray
            Corrected projection in um along z direction.
        """

        # Read the FEA file
        data = self.getMirrorData(FEAfileName, skiprows=1)

        # Z-gradient in um
        tzdz = data[:, 4]

        # r-gradient in um
        trdz = data[:, 5]

        # Get the temprature correction
        tempCorrInUm = M2TzGrad * tzdz + M2TrGrad * trdz

        return tempCorrInUm

    def getMirrorResInMmInZemax(self, gridFileName="M2_1um_grid.DAT",
                                numTerms=28, writeZcInMnToFilePath=None):
        """Get the residue of surface (mirror print along z-axis) in mm under
        the Zemax coordinate.

        This value is after the fitting with spherical Zernike polynomials
        (zk).
        
        Parameters
        ----------
        gridFileName : str, optional
            File name of bending mode data. (the default is "M2_1um_grid.DAT".)
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
        data = self.getMirrorData(gridFileName)

        # Get the x, y coordinate
        bx = data[:, 0]
        by = data[:, 1]

        # Transform the M2 coordinate to Zemax coordinate
        bxInZemax, byInZemax, surfInZemax = phosim2ZemaxCoorTrans(
                                                bx, by, self.getSurfAlongZ())

        # Get the mirror residue and zk in um
        RinM = self.getOuterRinM()

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

    def writeMirZkAndGridResInZemax(self, resFile="", surfaceGridN=200,
                                    gridFileName="M2_1um_grid.DAT",
                                    numTerms=28, writeZcInMnToFilePath=None):
        """Write the grid residue in mm of mirror surface after the fitting
        with Zk under the Zemax coordinate.

        Parameters
        ----------
        resFile : str, optional
            File path to save the grid surface residue map. (the default
            is "".)
        surfaceGridN : int, optional
            Surface grid number. (the default is 200.)
        gridFileName : str, optional
            File name of bending mode data. (the default is "M2_1um_grid.DAT".)
        numTerms : int, optional
            Number of Zernike terms to fit. (the default is 28.)
        writeZcInMnToFilePath : str, optional
            File path to write the fitted zk in mm. (the default is None.)

        Returns
        -------
        str
            Grid residue map related data.
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = \
            self.getMirrorResInMmInZemax(
                        gridFileName=gridFileName, numTerms=numTerms,
                        writeZcInMnToFilePath=writeZcInMnToFilePath)[0:3]

        # Change the unit from m to mm
        innerRinMm = self.getInnerRinM() * 1e3
        outerRinMm = self.getOuterRinM() * 1e3

        # Get the residue map used in Zemax
        # Content header: (NUM_X_PIXELS, NUM_Y_PIXELS, delta x, delta y)
        # Content: (z, dx, dy, dxdy)
        content = self._gridSampInMnInZemax(
                            resInMmInZemax, bxInMmInZemax, byInMmInZemax,
                            innerRinMm, outerRinMm, surfaceGridN, surfaceGridN,
                            resFile=resFile)

        return content

    def showMirResMap(self, gridFileName="M2_1um_grid.DAT", numTerms=28,
                      resFile=None, writeToResMapFilePath=None):
        """Show the mirror residue map.

        Parameters
        ----------
        gridFileName : str, optional
            File name of bending mode data. (the default is "M2_1um_grid.DAT".)
        numTerms : int, optional
            Number of Zernike terms to fit. (the default is 28.)
        resFile : str, optional
            File path of the grid surface residue map. (the default is None.)
        writeToResMapFilePath : str, optional
            File path to save the residue map. (the default is None.)
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = \
                    self.getMirrorResInMmInZemax(gridFileName=gridFileName,
                                                 numTerms=numTerms)[0:3]

        # Change the unit
        outerRinMm = self.getOuterRinM() * 1e3
        plotResMap(resInMmInZemax, bxInMmInZemax, byInMmInZemax,
                   outerRinMm, resFile=resFile,
                   writeToResMapFilePath=writeToResMapFilePath)


if __name__ == "__main__":
    pass
