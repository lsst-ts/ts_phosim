import os
import numpy as np
import unittest

from lsst.ts.phosim.PlotUtil import showFieldMap, plotFwhmOfIters
from lsst.ts.phosim.Utility import getModulePath


class TestPlotUtil(unittest.TestCase):
    """ Test the PlotUtil functions."""

    def setUp(self):

        modulePath = getModulePath()
        self.testData = os.path.join(modulePath, "tests", "testData")
        self.outFigFilePath = os.path.join(modulePath, "output", "img",
                                           "testFig.png")

    def tearDown(self):
        if (os.path.exists(self.outFigFilePath)):
            os.remove(self.outFigFilePath)

    def testShowFieldMap(self):

        fieldX = np.array([0, 0])
        fieldY = np.array([0, 1])

        showFieldMap(fieldX=fieldX, fieldY=fieldY,
                     saveToFilePath=self.outFigFilePath)
        self.assertTrue(os.path.exists(self.outFigFilePath))

    def testPlotFwhmOfIters(self):

        iterDataDir = os.path.join(self.testData, "iterData")
        pssnFiles = [os.path.join(iterDataDir, "iter%d" % num, "img", "PSSN.txt")
                     for num in range(5)]

        plotFwhmOfIters(pssnFiles, saveToFilePath=self.outFigFilePath)
        self.assertTrue(os.path.exists(self.outFigFilePath))


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
