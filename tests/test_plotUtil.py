import os
import numpy as np
import unittest

from lsst.ts.phosim.PlotUtil import showFieldMap
from lsst.ts.phosim.Utility import getModulePath


class TestPlotUtil(unittest.TestCase):
    """ Test the PlotUtil functions."""

    def testShowFieldMap(self):

        saveToFilePath = os.path.join(getModulePath(), "output", "img",
                                      "fieldMap.png")
        fieldX = np.array([0, 0])
        fieldY = np.array([0, 1])

        showFieldMap(fieldX=fieldX, fieldY=fieldY,
                     saveToFilePath=saveToFilePath)
        os.remove(saveToFilePath)


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
