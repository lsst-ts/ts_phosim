import os, unittest
import numpy as np

class SkySim(object):

    def __init__(self):
        """
        
        Initiate the SkySim object.
        """

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

    def resetSky(self):
        """
        
        Reset the sky information and delete all existed stars.
        """
        
        self.__init__()

    def addStarByFile(self, readFilePath, skiprows=1):
        """
        
        Add the star data by reading the file.
        
        Arguments:
            readFilePath {[star]} -- Star data file path.
        
        Keyword Arguments:
            skiprows {int} -- Skip the first "skiprows" lines. (default: {1})
        """
        
        # Read the data
        data = np.loadtxt(readFilePath, skiprows=skiprows)

        # Add the stars
        for star in data:
            self.addStarByRaDecInDeg(star[0], star[1], star[2], star[3])

    def addStarByChipPos(self):
        pass

    def addStarByFocalPlane(self):
        pass

class SkySimTest(unittest.TestCase):
    
    """
    Test functions in SkySim.
    """

    def setUp(self):

        self.skyFile = os.path.join("..", "data", "sky", "wfsStar.txt")

    def testFunc(self):

        # Instantiate the skySim object
        skySim = SkySim()

        skySim.setStarRaDecInDeg(np.array([0]), np.array([1]), np.array([2]), 
                                    np.array([3]))
        self.assertEqual(len(skySim.starId), 1)

        skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        self.assertEqual(len(skySim.starId), 2)

        try:
            skySim.addStarByRaDecInDeg(1, 2, 3, 4)
        except Exception as ValueError:
            pass

        skySim.resetSky()
        self.assertEqual(len(skySim.starId), 0)

        skySim.addStarByFile(self.skyFile)
        self.assertEqual(len(skySim.starId), 8)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
