import os
import numpy as np

from cwfs.Tool import ZernikeFit, ZernikeEval
# from lsst.cwfs.tools import ZernikeFit, ZernikeEval

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

    def gridSamp(self):
        pass

    def getLUTforce(self):
        pass

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

    def getMirrorResInZemax(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def getActForce(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def getPrintthz(self):
        raise NotImplementedError("Should have the child class implemented this.")

    def getTempCorr(self):
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

