import os
import argparse
import numpy as np

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import FilterType

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.Utility import getAoclcOutputPath, getConfigDir


def main(numStars, starSep, magFaintLim, magBrightLim,
         outputFilePath):

    """
    numStars: number of stars per field

    starSep: Minimum separation between stars

    magFaintLim: Faint Limit for magnitudes in catalog

    magBrightLim: Bright Limit for magnitudes in catalog

    outputFilePath: Filename for output catalog
    """

    # Survey parameters
    surveySettingFilePath = os.path.join(getConfigDir(),
                                         "surveySettings.yaml")
    surveySettings = ParamReader(filePath=surveySettingFilePath)
    filterType = FilterType.fromString(
        surveySettings.getSetting("filterType"))
    raInDeg = surveySettings.getSetting("raInDeg")
    decInDeg = surveySettings.getSetting("decInDeg")
    rotAngInDeg = surveySettings.getSetting("rotAngInDeg")

    ofcCalc = _prepareOfcCalc(filterType, rotAngInDeg)
    skySim = SkySim()
    metr = OpdMetrology()
    metr.setDefaultComcamGQ()

    starMagList = np.linspace(magBrightLim, magFaintLim, numStars)
    skySim = _addStarsInField(skySim, metr, numStars,
                              starSep, starMagList)
    skySim.exportSkyToFile(outputFilePath)


def _prepareOfcCalc(filterType, rotAngInDeg):

    ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
    ofcCalc.setFilter(filterType)
    ofcCalc.setRotAng(rotAngInDeg)
    ofcCalc.setGainByPSSN()

    return ofcCalc


def _addStarsInField(skySim, opdMetr, numStars, starSep, starMagList):

    starId = 0
    raInDegList, declInDegList = opdMetr.getFieldXY()

    raOffset = np.zeros(numStars)
    decMin = -1 * starSep * (float(numStars-1) / 2)
    decMax = -1 * decMin
    decOffset = np.linspace(decMin, decMax, numStars)

    for raInDeg, declInDeg in zip(raInDegList, declInDegList):
        # It is noted that the field position might be < 0. But it is not the
        # same case for ra (0 <= ra <= 360).
        if (raInDeg < 0):
            raInDeg += 360.0

        for num in range(numStars):
            skySim.addStarByRaDecInDeg(starId, raInDeg+raOffset[num],
                                       declInDeg+decOffset[num], starMagList[num])
            starId += 1

    return skySim


if __name__ == "__main__":

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Create a Star Catalog for comcamCloseLoop simulation.")
    parser.add_argument("--numStars", type=int, default=1,
                        help="number of stars per field (default: 1)")
    parser.add_argument("--starSep", type=float, default=0.05,
                        help="minimum separation between stars in degrees (default: 0.05)")
    parser.add_argument("--magFaintLim", type=float, default=15.,
                        help="faint limit of stars in catalog (default: 15)")
    parser.add_argument("--magBrightLim", type=float, default=15.,
                        help="bright limit of stars in catalog (default: 15)")
    parser.add_argument("--output", type=str, default="",
                        help="output catalog filepath")
    args = parser.parse_args()

    if (args.output == ""):
        outputDir = getAoclcOutputPath()
        outputFilePath = os.path.join(outputDir, "starCat.txt")
    else:
        outputFilePath = args.output

    main(args.numStars, args.starSep, args.magFaintLim, args.magBrightLim,
         outputFilePath)
