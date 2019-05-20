import os
import unittest

from lsst.ts.phosim.Utility import opt2ZemaxCoorTrans, zemax2optCoorTrans, \
    mapSurfNameToEnum, SurfaceType, getPhoSimPath, sortOpdFileList, \
    getAoclcOutputPath, getModulePath


class TestUtility(unittest.TestCase):
    """Test the Utility functions."""

    def testPhosim2ZemaxCoorTrans(self):

        xInPhosim, yInPhosim, zInPhosim = 1, 2, 3
        xInZemax, yInZemax, zInZemax = opt2ZemaxCoorTrans(
            xInPhosim, yInPhosim, zInPhosim)
        self.assertEqual((xInZemax, yInZemax, zInZemax),
                         (-xInPhosim, yInPhosim, -zInPhosim))

    def testZemax2phosimCoorTrans(self):

        xInZemax, yInZemax, zInZemax = 1, 2, 3
        xInPhosim, yInPhosim, zInPhosim = zemax2optCoorTrans(
            xInZemax, yInZemax, zInZemax)
        self.assertEqual((xInPhosim, yInPhosim, zInPhosim),
                         (-xInZemax, yInZemax, -zInZemax))

    def testMapSurfNameToEnum(self):

        surfName = "L2F"
        surfaceType = mapSurfNameToEnum(surfName)
        self.assertEqual(surfaceType, SurfaceType.L2F)
        self.assertRaises(ValueError, mapSurfNameToEnum, "L123F")

    def testGetPhoSimPathNotExist(self):

        self.assertRaises(RuntimeError, getPhoSimPath, "WRONGPATH")

    def testGetPhoSimPath(self):

        PHOSIMPATH = "/path/to/phosim"
        os.environ["PHOSIMPATH"] = PHOSIMPATH

        phosimPath = getPhoSimPath()
        self.assertEqual(phosimPath, PHOSIMPATH)

        os.environ.pop("PHOSIMPATH")

    def testGetAoclcOutputPathNotAssigned(self):

        aoclcOutputPath = getAoclcOutputPath()
        self.assertEqual(aoclcOutputPath,
                         os.path.join(getModulePath(), "output"))

    def testGetAoclcOutputPath(self):

        AOCLCOUTPUTPATH = "/path/to/aoclc/output"
        os.environ["AOCLCOUTPUTPATH"] = AOCLCOUTPUTPATH

        aoclcOutputPath = getAoclcOutputPath()
        self.assertEqual(aoclcOutputPath, AOCLCOUTPUTPATH)

        os.environ.pop("AOCLCOUTPUTPATH")

    def testSortOpdFileList(self):

        fileDir = "/fileDir"
        opdFileNameList = ["opd_100_3.fits.gz", "opd_100_1.fits.gz",
                           "opd_100_13.fits.gz", "opd_100_0.fits.gz",
                           "opd_100_10.fits.gz"]
        opdFileList = [os.path.join(fileDir, fileName)
                       for fileName in opdFileNameList]

        sortedOpdFileList = sortOpdFileList(opdFileList)

        sortedIdxList = [3, 1, 0, 4, 2]
        ansOpdFileList = [opdFileList[idx] for idx in sortedIdxList]
        self.assertListEqual(sortedOpdFileList, ansOpdFileList)

    def testSortOpdFileListWithUnmatchedName(self):

        self.assertRaises(ValueError, sortOpdFileList, ["wrong_name"])


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
