import os
import numpy as np

from lsst.ts.phosim.telescope.MirrorSim import MirrorSim

from lsst.ts.phosim.Utility import opt2ZemaxCoorTrans
from lsst.ts.phosim.PlotUtil import plotResMap
from lsst.ts.phosim.Utility import getConfigDir

from lsst.ts.wep.ParamReader import ParamReader


class M2Sim(MirrorSim):

    def __init__(self):
        """Initiate the M2 simulator class."""

        # M2 setting file
        configDir = os.path.join(getConfigDir(), "M2")
        settingFilePath = os.path.join(configDir, "m2Setting.yaml")
        self._m2SettingFile = ParamReader(filePath=settingFilePath)

        # Inner and outer radius of M2 mirror in m
        radiusInner = self._m2SettingFile.getSetting("radiusInner")
        radiusOuter = self._m2SettingFile.getSetting("radiusOuter")

        super(M2Sim, self).__init__(radiusInner, radiusOuter, configDir)

        # Mirror surface bending mode grid file
        self._gridFile = ParamReader()

        # Mirror FEA model with gradient temperature data
        self._feaFile = ParamReader()

        self._config("M2_1um_force.yaml", "", "M2_1um_grid.yaml",
                     "M2_GT_FEA.yaml")

    def _config(self, actForceFileName, lutFileName, gridFileName, feaFileName):
        """Do the configuration.

        LUT: Look-up table.
        FEA: Finite element analysis.

        Parameters
        ----------
        actForceFileName : str
            Actuator force file name.
        lutFileName : str
            LUT file name.
        gridFileName : str
            File name of bending mode data.
        feaFileName : str
            FEA model data file name.
        """

        numTerms = self._m2SettingFile.getSetting("numTerms")

        super(M2Sim, self).config(numTerms=numTerms,
                                  actForceFileName=actForceFileName,
                                  lutFileName=lutFileName)

        mirrorDataDir = self.getMirrorDataDir()

        gridFilePath = os.path.join(mirrorDataDir, gridFileName)
        self._gridFile.setFilePath(gridFilePath)

        feaFilePath = os.path.join(mirrorDataDir, feaFileName)
        self._feaFile.setFilePath(feaFilePath)

    def getPrintthz(self, zAngleInRadian, preCompElevInRadian=0):
        """Get the mirror print in um along z direction in specific zenith
        angle.

        FEA: Finite element analysis.

        Parameters
        ----------
        zAngleInRadian : float
            Zenith angle in radian.
        preCompElevInRadian : float, optional
            Pre-compensation elevation angle in radian. (the default is 0.)

        Returns
        ------
        numpy.ndarray
            Corrected projection in um along z direction.
        """

        # Read the FEA file
        data = self._feaFile.getMatContent()

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

    def getTempCorr(self, m2TzGrad, m2TrGrad):
        """Get the mirror print correction along z direction for certain
        temperature gradient.

        FEA: Finite element analysis.

        Parameters
        ----------
        m2TzGrad : float
            Temperature gradient along z direction in degree C (+/-2sigma
            spans 1C).
        m2TrGrad : float
            Temperature gradient along r direction in degree C (+/-2sigma
            spans 1C).

        Returns
        -------
        numpy.ndarray
            Corrected projection in um along z direction.
        """

        # Read the FEA file
        data = self._feaFile.getMatContent()

        # Z-gradient in um
        tzdz = data[:, 4]

        # r-gradient in um
        trdz = data[:, 5]

        # Get the temprature correction
        tempCorrInUm = m2TzGrad * tzdz + m2TrGrad * trdz

        return tempCorrInUm

    def getMirrorResInMmInZemax(self, writeZcInMnToFilePath=None):
        """Get the residue of surface (mirror print along z-axis) in mm under
        the Zemax coordinate.

        This value is after the fitting with spherical Zernike polynomials
        (zk).

        Parameters
        ----------
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
        data = self._gridFile.getMatContent()

        # Get the x, y coordinate
        bx = data[:, 0]
        by = data[:, 1]

        # Transform the M2 coordinate to Zemax coordinate
        bxInZemax, byInZemax, surfInZemax = opt2ZemaxCoorTrans(
            bx, by, self.getSurfAlongZ())

        # Get the mirror residue and zk in um
        RinM = self.getOuterRinM()
        resInUmInZemax, zcInUmInZemax = self._getMirrorResInNormalizedCoor(
            surfInZemax, bxInZemax/RinM, byInZemax/RinM)

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
                                    writeZcInMnToFilePath=None):
        """Write the grid residue in mm of mirror surface after the fitting
        with Zk under the Zemax coordinate.

        Parameters
        ----------
        resFile : str, optional
            File path to save the grid surface residue map. (the default
            is "".)
        surfaceGridN : int, optional
            Surface grid number. (the default is 200.)
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
                writeZcInMnToFilePath=writeZcInMnToFilePath)[0:3]

        # Change the unit from m to mm
        innerRinMm = self.getInnerRinM() * 1e3
        outerRinMm = self.getOuterRinM() * 1e3

        # Get the residue map used in Zemax
        # Content header: (NUM_X_PIXELS, NUM_Y_PIXELS, delta x, delta y)
        # Content: (z, dx, dy, dxdy)
        content = self._gridSampInMnInZemax(
            resInMmInZemax, bxInMmInZemax, byInMmInZemax, innerRinMm,
            outerRinMm, surfaceGridN, surfaceGridN, resFile=resFile)

        return content

    def showMirResMap(self, resFile, writeToResMapFilePath=None):
        """Show the mirror residue map.

        Parameters
        ----------
        resFile : str
            File path of the grid surface residue map.
        writeToResMapFilePath : str, optional
            File path to save the residue map. (the default is None.)
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = \
            self.getMirrorResInMmInZemax()[0:3]

        # Change the unit
        outerRinMm = self.getOuterRinM() * 1e3
        plotResMap(resInMmInZemax, bxInMmInZemax, byInMmInZemax,
                   outerRinMm, resFile=resFile,
                   writeToResMapFilePath=writeToResMapFilePath)


if __name__ == "__main__":
    pass
