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


class createPhosimCatalog():

    def createPhosimCatalog(self, numStars, starSep, magList,
                            raOffset, decOffset, outputFilePath):

        """
        numStars: number of stars per field

        starSep: Minimum separation between stars

        magList: Star magnitudes in catalog

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

        ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg)
        skySim = SkySim()
        metr = OpdMetrology()
        metr.setDefaultComcamGQ()

        skySim = self._addStarsInField(skySim, metr, numStars,
                                       starSep, raOffset, decOffset,
                                       magList)
        skySim.exportSkyToFile(outputFilePath)


    def _prepareOfcCalc(self, filterType, rotAngInDeg):

        ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
        ofcCalc.setFilter(filterType)
        ofcCalc.setRotAng(rotAngInDeg)
        ofcCalc.setGainByPSSN()

        return ofcCalc


    def _addStarsInField(self, skySim, opdMetr, numStars, starSep,
                         raOffset, decOffset, starMagList):

        starId = 0
        raInDegList, declInDegList = opdMetr.getFieldXY()

        raShift = np.zeros(numStars)
        decMin = -1 * starSep * (float(numStars-1) / 2)
        decMax = -1 * decMin
        decShift = np.linspace(decMin, decMax, numStars)

        raShift += raOffset
        decShift += decOffset

        for raInDeg, declInDeg in zip(raInDegList, declInDegList):
            # It is noted that the field position might be < 0. But it is not the
            # same case for ra (0 <= ra <= 360).
            # if (raInDeg < 0):
            #     raInDeg += 360.0

            for num in range(numStars):
                raPerturbed = raInDeg + raShift[num]
                decPerturbed = declInDeg + decShift[num]

                if raPerturbed < 0:
                    raPerturbed += 360.0

                skySim.addStarByRaDecInDeg(starId, raPerturbed,
                                           decPerturbed, starMagList[num])
                starId += 1

        return skySim


if __name__ == "__main__":

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Create a Star Catalog for comcamCloseLoop simulation.")
    parser.add_argument("--numStars", type=int, default=1,
                        help="number of stars per field (default: 1)")
    parser.add_argument("--raOffset", type=float, default=0.0,
                        help="Offset from center of chip in RA in arcsec (default: 0.)")
    parser.add_argument("--decOffset", type=float, default=0.0,
                        help="Offset from center of chip in dec in arcsec (default: 0.)")
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

    starMagList = np.linspace(args.magBrightLim, args.magFaintLim, args.numStars)

    createCat = createPhosimCatalog()
    createCat.createPhosimCatalog(args.numStars, args.starSep,
                                  args.magFaintLim, args.magBrightLim,
                                  args.raOffset, args.decOffset,
                                  outputFilePath)
