import os
import numpy as np
from scipy.interpolate import Rbf
import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lsst.ts.wep.cwfs.Tool import ZernikeFit, ZernikeEval


class MirrorSim(object):
    
    def __init__(self, innerRinM, outerRinM, mirrorDataDir=""):
        """Initiate the MirrorSim object.

        Parameters
        ----------
        innerRinM : float
            Mirror inner radius in m.
        outerRinM : float
            Mirror outer radius in m.
        mirrorDataDir : str, optional
            Mirror data directory. (the default is "".)
        """

        self.RiInM = innerRinM
        self.RinM = outerRinM
        self.mirrorDataDir = mirrorDataDir

        self._surf = np.array([])

    def setMirrorDataDir(self, mirrorDataDir):
        """Set the directory of mirror data.

        Parameters
        ----------
        mirrorDataDir : str
            Directory to mirror data.
        """

        self.mirrorDataDir = mirrorDataDir

    def getMirrorData(self, dataFileName, skiprows=0):
        """Get the mirror data.

        Parameters
        ----------
        dataFileName : str
            Data file name.
        skiprows : int, optional
            Skip the first 'skiprows' lines (the default is 0.)
        
        Returns
        -------
        numpy.ndarray
            Mirror data.
        """

        filePath = os.path.join(self.mirrorDataDir, dataFileName)
        data = np.loadtxt(filePath, skiprows=skiprows)

        return data

    def setSurfAlongZ(self, surfAlongZinUm):
        """Set the mirror surface along the z direction in um.

        Parameters
        ----------
        surfAlongZinUm : numpy.ndarray
            Mirror surface along the z direction in um.
        """

        self._surf = np.array(surfAlongZinUm, dtype=float)

    def getSurfAlongZ(self):
        """Get the mirror surface along the z direction in um.

        Returns
        -------
        numpy.ndarray
            Mirror surface along the z direction in um.
        """

        return self._surf

    def getLUTforce(self, zangleInDeg, LUTfileName):
        """Get the actuator force of mirror based on LUT.
        
        LUT: Look-up table.
        
        Parameters
        ----------
        zangleInDeg : float
            Zenith angle in degree.
        LUTfileName : str
            LUT file name.
        
        Returns
        -------
        numpy.ndarray
            Actuator forces in specific zenith angle.
        
        Raises
        ------
        ValueError
            The degee order in LUT is incorrect.
        """

        # Read the LUT file
        lut = self.getMirrorData(LUTfileName)

        # Get the step. The values of LUT are listed in every step size.
        # The degree range is 0 - 90 degree.
        # The file in the simulation is every 1 degree. The formal one should 
        # be every 5 degree.
        ruler = lut[0, :]
        stepList = np.diff(ruler)
        if np.any(stepList <= 0):
            raise ValueError("The degee order in LUT is incorrect.")

        # If the specific zenith angle is larger than the listed angle range,
        # use the biggest listed zenith angle data instead.
        if (zangleInDeg >= ruler.max()):
            lutForce = lut[1:, -1]

        # If the specific zenith angle is smaller than the listed angle range,
        # use the smallest listed zenith angle data instead.
        elif (zangleInDeg <= ruler.min()):
            lutForce = lut[1:, 0]

        # If the specific zenith angle is in the listed angle range,
        # do the linear fit to get the data.
        else:
            # Find the boundary indexes for the specific zenith angle
            p1 = np.where(ruler<=zangleInDeg)[0][-1]
            p2 = p1+1

            # Do the linear approximation
            w2 = (zangleInDeg-ruler[p1])/stepList[p1]
            w1 = 1-w2

            lutForce = w1*lut[1:, p1] + w2*lut[1:, p2]

        return lutForce

    def __gridSampInMnInZemax(self, zfInMm, xfInMm, yfInMm, innerRinMm, outerRinMm, nx, ny, resFile=None):
        """
        
        Get the grid residue map used in Zemax.
        
        Arguments:
            zfInMm {[ndarray]} -- Surface map in mm.
            xfInMm {[ndarray]} -- x position in mm.
            yfInMm {[ndarray]} -- y position in mm.
            innerRinMm {[float]} -- Inner radius in mm.
            outerRinMm {[float]} -- Outer radius in mm.
            nx {[int]} -- Number of pixel along x-axis of surface residue map. It is noted that 
                          the real pixel number is nx + 4.
            ny {[int]} -- Number of pixel along y-axis of surface residue map. It is noted that 
                          the real pixel number is ny + 4. 

        Keyword Arguments:
            resFile {[str]} -- File path to write the surface residue map. (default: {None})

        Returns:
            [str] -- Grid residue map related data.
        """

        # Radial basis function approximation/interpolation of surface
        Ff = Rbf(xfInMm, yfInMm, zfInMm)

        # Number of grid points on x-, y-axis. 
        # Alway extend 2 points on each side
        # Do not want to cover the edge? change 4->2 on both lines
        NUM_X_PIXELS = nx+4
        NUM_Y_PIXELS = ny+4

        # This is spatial extension factor, which is calculated by the slope at edge
        extFx = (NUM_X_PIXELS-1) / (nx-1)
        extFy = (NUM_Y_PIXELS-1) / (ny-1)
        extFr = np.sqrt(extFx * extFy)

        # Delta x and y 
        delx = outerRinMm*2*extFx / (NUM_X_PIXELS-1)
        dely = outerRinMm*2*extFy / (NUM_Y_PIXELS-1)

        # Minimum x and y
        minx = -0.5*(NUM_X_PIXELS-1)*delx
        miny = -0.5*(NUM_Y_PIXELS-1)*dely

        # Calculate the epsilon
        epsilon = 1e-4*min(delx, dely)
        
        # Write four numbers for the header line
        content = "%d %d %.9E %.9E\n" % (NUM_X_PIXELS, NUM_Y_PIXELS, delx, dely)

        #  Write the rows and columns
        for jj in range(1, NUM_X_PIXELS + 1):
            for ii in range(1, NUM_Y_PIXELS + 1):

                # x and y positions
                x =  minx + (ii - 1) * delx
                y =  miny + (jj - 1) * dely
                
                # Invert top to bottom, because Zemax reads (-x,-y) first
                y = -y

                # Calculate the radius
                r = np.sqrt(x**2 + y**2)

                # Set the value as zero when the radius is not between the inner and outer radius.
                if (r < innerRinMm/extFr) or (r > outerRinMm*extFr):
                    
                    z = 0
                    dx = 0
                    dy = 0
                    dxdy = 0

                # Get the value by the fitting
                else:

                    # Get the z
                    z = Ff(x, y)
                    
                    # Compute the dx
                    tem1 = Ff((x+epsilon), y)
                    tem2 = Ff((x-epsilon), y)
                    dx = (tem1 - tem2)/(2.0*epsilon)

                    # Compute the dy
                    tem1 = Ff(x, (y+epsilon))
                    tem2 = Ff(x, (y-epsilon))
                    dy = (tem1 - tem2)/(2.0*epsilon)

                    # Compute the dxdy
                    tem1 = Ff((x+epsilon), (y+epsilon))
                    tem2 = Ff((x-epsilon), (y+epsilon))
                    tem3 = (tem1 - tem2)/(2.0*epsilon)
                    
                    tem1 = Ff((x+epsilon), (y-epsilon))
                    tem2 = Ff((x-epsilon), (y-epsilon))
                    tem4 = (tem1 - tem2)/(2.0*epsilon)
                    
                    dxdy = (tem3 - tem4)/(2.0*epsilon)

                content += "%.9E %.9E %.9E %.9E\n" % (z, dx, dy, dxdy)

        # Write the surface residue data into the file
        if (resFile is not None):
            outid = open(resFile, "w");
            outid.write(content)
            outid.close()

        return content

    def __showResMap(self, zfInMm, xfInMm, yfInMm, outerRinMm, resFile=None, writeToResMapFilePath=None):
        """
        
        Show the mirror residue map.
        
        Arguments:
            zfInMm {[ndarray]} -- Surface map in mm.
            xfInMm {[ndarray]} -- x position in mm.
            yfInMm {[ndarray]} -- y position in mm.
            outerRinMm {[float]} -- Outer radius in mm.
        
        Keyword Arguments:
            resFile {[str]} -- File path of the grid surface residue map. (default: {None})
            writeToResMapFilePath {[str]} -- File path to save the residue map. (default: {None})
        """

        # Plot the figure
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))

        # The input data to gridSamp.m is in mm (zemax default)
        sc = ax[1].scatter(xfInMm, yfInMm, s=25, c=zfInMm*1e6, marker=".", edgecolor="none")
        ax[1].axis("equal")
        ax[1].set_title("Surface map on FEA grid (nm)")
        ax[1].set_xlim([-outerRinMm, outerRinMm])
        ax[1].set_ylim([-outerRinMm, outerRinMm])
        ax[1].set_xlabel("x (mm)")

        if (resFile is not None):

            # Get the data
            data = np.loadtxt(resFile)
            NUM_X_PIXELS, NUM_Y_PIXELS, delxInMm, delyInMm = data[0,:]
            NUM_X_PIXELS = int(NUM_X_PIXELS)
            NUM_Y_PIXELS = int(NUM_Y_PIXELS)

            # Get the zp data
            zpTemp = data[1:, 0]
            zp = np.zeros((NUM_X_PIXELS, NUM_Y_PIXELS))
            for jj in range(1, NUM_X_PIXELS + 1):
                for ii in range(1, NUM_Y_PIXELS + 1):
                    zp[NUM_X_PIXELS+1-jj-1, ii-1] = zpTemp[(jj-1)*NUM_X_PIXELS + (ii-1)]

            # Minimum x and y
            minx = -0.5*(NUM_X_PIXELS-1)*delxInMm
            miny = -0.5*(NUM_Y_PIXELS-1)*delyInMm

            xx = np.linspace(minx, -minx, NUM_X_PIXELS)
            yy = np.linspace(miny, -miny, NUM_Y_PIXELS)
            xp, yp = np.meshgrid(xx, yy)

            xp = xp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS, 1))
            yp = yp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS, 1))
            zp = zp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS, 1))

            sc = ax[0].scatter(xp, yp, s=25, c=zp*1e6, marker=".", edgecolor="none")

        ax[0].axis("equal")
        ax[0].set_title("grid input to ZEMAX (nm)")
        ax[0].set_xlim([-outerRinMm, outerRinMm])
        ax[0].set_ylim([-outerRinMm, outerRinMm])
        ax[0].set_xlabel("x (mm)")
        ax[0].set_ylabel("y (mm)")

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        fig.colorbar(sc, cax=cbar_ax)

        if (writeToResMapFilePath is not None):
            plt.savefig(writeToResMapFilePath)
            plt.close()
        else:
            plt.show()

    def __getMirrorResInNormalizedCoor(self, surf, x, y, numTerms):
        """
        
        Get the residue of surface (mirror print along z-axis) after the fitting with Zk in the 
        normalized x, y coordinate.
        
        Arguments:
            x {[ndarray]} -- Normalized x coordinate.
            y {[ndarray]} -- Normalized y coordinate.
            numTerms {[int]} -- Number of Zernike terms to fit.
                
        Returns:
            [ndarray] -- Surface residue after the fitting.
            [ndarray] -- Fitted Zernike polynomials.
        """

        # Get the surface change along the z-axis in the basis of Zk
        # It is noticed that the x and y coordinates are normalized for the fitting 
        zc = ZernikeFit(surf, x, y, numTerms)

        # Residue of fitting
        res = surf - ZernikeEval(zc, x, y)

        return res, zc

    def getActForce(self):
        raise NotImplementedError("Child class should implemented this.")

    def getPrintthz(self):
        raise NotImplementedError("Child class should implemented this.")

    def getTempCorr(self):
        raise NotImplementedError("Child class should implemented this.")

    def getMirrorResInMmInZemax(self):
        raise NotImplementedError("Child class should implemented this.")

    def writeMirZkAndGridResInZemax(self):
        raise NotImplementedError("Child class should implemented this.")

    def showMirResMap(self):
        raise NotImplementedError("Child class should implemented this.")


if __name__ == "__main__":

    # Do the unit test
    unittest.main()
