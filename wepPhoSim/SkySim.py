import os
import numpy as np

class SkySim(object):

    def __init__(self):

        self.starId = np.array([])
        self.ra = np.array([])
        self.decl = np.array([])
        self.mag = np.array([])

    def setStarRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """
        
        Set the star information by (ra, dec) in degrees.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
            raInDeg {[list/ ndarray]} -- Star Ra in degree.
            declInDeg {[list/ ndarray]} -- Star Decl in degree.
            mag {[list/ ndarray]} -- Star magnitude.
        """

        self.__setUniqStarId(starId)
        self.ra = raInDeg
        self.decl = declInDeg
        self.mag = mag

    def __setUniqStarId(self, starId):
        """
        
        Set the star unique IDs.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
        
        Raises:
            ValueError -- Star IDs are not unique.
        """

        # Collect all star id
        allStarId = np.append(self.starId, starId).astype("int")

        # Check all IDs are unique all not
        if (len(allStarId) != len(np.unique(allStarId))):
            raise ValueError("Star IDs are not unique.")

        # Set the star ID
        self.starId = allStarId

    def addStarByRaDecInDeg(self, starId, raInDeg, declInDeg, mag):
        """
        
        Add the star information by (ra, dec) in degrees.
        
        Arguments:
            starId {[list/ ndarray]} -- Star ID.
            raInDeg {[list/ ndarray]} -- Star Ra in degree.
            declInDeg {[list/ ndarray]} -- Star Decl in degree.
            mag {[list/ ndarray]} -- Star magnitude.
        """

        self.__setUniqStarId(starId)
        self.ra = np.append(self.ra, raInDeg)
        self.decl = np.append(self.decl, declInDeg)
        self.mag = np.append(self.mag, mag)

    def addStarByFile(self):
        pass

    def addStarByChipPos(self):
        pass

    def addStarByFocalPlane(self):
        pass

if __name__ == "__main__":

    skySim = SkySim()

    starId = 0
    raInDeg = 0.1
    declInDeg = 0.1
    skySim.addStarByRaDecInDeg(starId, raInDeg, declInDeg)