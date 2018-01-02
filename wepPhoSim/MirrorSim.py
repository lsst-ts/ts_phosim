import os
import numpy as np
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt

# from cwfs.Tool import ZernikeFit, ZernikeEval
from lsst.cwfs.tools import ZernikeFit, ZernikeEval

class MirrorSim(object):
    
    def __init__(self, innerRinM, outerRinM, surf=None, mirrorDataDir=None):

        # Check the raidus
        if (innerRinM >= outerRinM):
            raise ValueError("Inner raidus should be smaller than outer raidus.")
        
        self.RiInM = innerRinM
        self.RinM = outerRinM
        self.surf = surf

        self.mirrorDataDir = mirrorDataDir

    def setMirrorDataDir(self, mirrorDataDir):
        """
        
        Set the directory of mirror data.
        
        Arguments:
            mirrorDataDir {[str]} -- Directory to mirror data.
        """

        self.mirrorDataDir = mirrorDataDir

    def getMirrorData(self, dataFileName, skiprows=0):
        """
        
        Get the mirror data.
        
        Arguments:
            dataFileName {[str]} -- Data file name.
        
        Keyword Arguments:
            skiprows {int} -- Skip the first "skiprows" lines. (default: {0})
        
        Returns:
            [ndarray] -- Mirror data.
        """

        data = np.loadtxt(os.path.join(self.mirrorDataDir, dataFileName), skiprows=skiprows)

        return data

    def setSurfAlongZ(self, surfAlongZ):
        """
        
        Set the mirror surface along the z direction.
        
        Arguments:
            surfAlongZ {[ndarray]} -- Mirror surface along the z direction.
        """

        self.surf = surfAlongZ

    def writeSurfRes(self):
        pass

    def getLUTforce(self):
        pass

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

    def __showResMap(self, zp, zfInMm, xfInMm, yfInMm, outerRinMm, delxInMm, delyInMm, NUM_X_PIXELS, NUM_Y_PIXELS, resFile):

        # Parameters

        # Minimum x and y
        minx = -0.5*(NUM_X_PIXELS-1)*delxInMm
        miny = -0.5*(NUM_Y_PIXELS-1)*delyInMm

        # Plot the figure
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))

        # The input data to gridSamp.m is in mm (zemax default)
        sc = ax[1].scatter(xf, yf, s=25, c=zf*1e6, marker=".", edgecolor="none")
        ax[1].axis("equal")
        ax[1].set_title("Surface map on FEA grid (nm)")
        ax[1].set_xlim([-outerR, outerR])
        ax[1].set_ylim([-outerR, outerR])
        ax[1].set_xlabel("x (mm)")

        xx = np.arange(minx, -minx + delx, delx)
        yy = np.arange(miny, -miny + dely, dely)
        xp, yp = np.meshgrid(xx, yy)
        xp = xp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS,1))
        yp = yp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS,1))
        zp = zp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS,1))
        sc = ax[0].scatter(xp, yp, s=25, c=zp*1e6, marker=".", edgecolor="none")

        ax[0].axis("equal")
        ax[0].set_title("grid input to ZEMAX (nm)")
        ax[0].set_xlim([-outerR, outerR])
        ax[0].set_ylim([-outerR, outerR])
        ax[0].set_xlabel("x (mm)")
        ax[0].set_ylabel("y (mm)")

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        fig.colorbar(sc, cax=cbar_ax)

        plt.savefig(resFile.replace(".txt", ".png"))

    def __getMirrorResInNormalizedCoor(self, surf, x, y, numTerms, writeZcToFilePath=None):
        """
        
        Get the residue of surface (mirror print along z-axis) after the fitting with Zk in the 
        normalized x, y coordinate.
        
        Arguments:
            x {[ndarray]} -- Normalized x coordinate.
            y {[ndarray]} -- Normalized y coordinate.
            numTerms {[int]} -- Number of Zernike terms to fit.
        
        Keyword Arguments:
            writeZcToFilePath {[str]} -- File path to save the fitted Zk. (default: {None})
        
        Returns:
            [ndarray] -- Surface residue after the fitting.
            [ndarray] -- Fitted Zernike polynomials.
        """

        # Get the surface change along the z-axis in the basis of Zk
        # It is noticed that the x and y coordinates are normalized for the fitting 
        zc = ZernikeFit(surf, x, y, numTerms)

        # Residue of fitting
        res = surf - ZernikeEval(zc, x, y)

        # Save the file of fitted Zk
        if (writeZcToFilePath is not None):
            np.savetxt(writeZcToFilePath, zc)

        return res, zc

    def getActForce(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def getPrintthz(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def getTempCorr(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def getMirrorResInZemax(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def writeMirZkAndGridResInZemax(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def showMirResMap(self):
        raise NotImplementedError("Should have the child class implemented this.")

if __name__ == "__main__":
    
    # Inner raidus
    innerRinM = 0.9

    # Outer raidus
    outerRinM = 1.710

    # M2 data directory
    mirrorDataDir = "../data/M2"
    actForceFileName = "M2_1um_force.DAT"

    # Instantiate mirror
    M2 = MirrorSim(innerRinM, outerRinM)
    M2.setMirrorDataDir(mirrorDataDir)

