import os
import numpy as np

class SkySim(object):

    def __init__(self):

        self.starId = np.array([])
        self.ra = np.array([])
        self.decl = np.array([])
        self.mag = np.array([])

    def setStarRaDec(self, starId, raInDeg, declInDeg, mag):

        self.__setUniqStarId(starId)
        self.ra = raInDeg
        self.decl = declInDeg
        self.mag = mag

    def addStarByRaDecInDeg(self, starId, raInDeg, declInDeg, mag):

        self.__setUniqStarId(starId)
        self.ra = np.append(self.ra, raInDeg)
        self.decl = np.append(self.decl, declInDeg)
        self.mag = np.append(self.mag, mag)

    def __setUniqStarId(self, starId):

        # Collect all star id
        allStarId = np.append(self.starId, starId).astype("int")

        # Check all IDs are unique all not
        if (len(allStarId) != len(np.unique(allStarId))):
            raise ValueError("Star IDs are not unique.")

        # Set the star ID
        self.starId = allStarId

if __name__ == "__main__":

    skySim = SkySim()

    starId = 0
    raInDeg = 0.1
    declInDeg = 0.1
    skySim.addStarByRaDecInDeg(starId, raInDeg, declInDeg)