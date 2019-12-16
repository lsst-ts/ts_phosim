import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from lsst.ts.phosim.Utility import getAoclcOutputPath


class calcMetrics():

    def loadZernikeData(self, opdZkFilePath, wfsZkFilePath):

        opdZkData = np.loadtxt(opdZkFilePath)
        wfsZkData = np.loadtxt(wfsZkFilePath)

        return opdZkData, wfsZkData


    def calcSSR(self, opdZkData, wfsZkData):

        zerDiffSqByField = np.sum((wfsZkData - opdZkData)**2, axis=1)
        zerDiffSqTotal = np.sum(zerDiffSqByField)

        return zerDiffSqByField, zerDiffSqTotal

    def plotZernikeDiff(self, opdZkData, wfsZkData, testLabel, saveToFilePath=None, dpi=None):
        """Plot the Difference between Wavefront Estimation Zernickes and OPD.

        Parameters
        ----------
        opdZkData:
            Zernike polynomials from OPD. Each row are the polynomials
            for a different field.
        wfsZkData:
            Zernike polynomials from Wavefront Estimation.
            Each row are the polynomials for a different field.
        saveToFilePath : str, optional
            File path to save the figure. If None, the figure will be showed. (the
            default is None.)
        dpi : int, optional
            The resolution in dots per inch. (the default is None.)
        """

        # Collect the data. Each row is the set Zernicke polynomials for a field.

        if len(opdZkData) != len(wfsZkData):
            raise IOError("Zernike Files do not have same number of rows.")

        zerDiff = wfsZkData - opdZkData
        numOfFields = len(opdZkData)

        # Plot the figure
        plt.figure()
        for row_num in range(numOfFields):
            plt.plot(np.arange(4, 23), zerDiff[row_num], "x-", label="Field %i" % (row_num+1))
        plt.xlabel("Zernike Polynomial Number")
        plt.ylabel("Difference from OPD")
        plt.legend()

        plt.title(testLabel)

        plt.savefig(saveToFilePath, dpi=dpi)

if __name__ == "__main__":

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run analysis comparing Zernike polynomials from WEP to OPD.")
    parser.add_argument("--input", type=str, default="",
                        help="input directory")
    parser.add_argument("--testLabel", type=str, default="")

    args = parser.parse_args()

    if (args.input == ""):
        inputDir = getAoclcOutputPath()
    else:
        inputDir = args.input

    calcMetrics(inputDir, args.testLabel)