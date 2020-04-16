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


class createBlendTestCatalog():

    def createPhosimCatalog(self, testStarInputFile, outputFilePath,
                            numFields=9):

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

        skySim = self._addStarsInField(skySim, metr, testStarInputFile)
        skySim.exportSkyToFile(outputFilePath)


    def _prepareOfcCalc(self, filterType, rotAngInDeg):

        ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
        ofcCalc.setFilter(filterType)
        ofcCalc.setRotAng(rotAngInDeg)
        ofcCalc.setGainByPSSN()

        return ofcCalc

    def _addStarsInField(self, skySim, opdMetr, testStarInputFile):

        raInDegList, declInDegList = opdMetr.getFieldXY()

        rand_state = np.random.RandomState(seed=413)

        blend_input = np.genfromtxt(testStarInputFile, names=True)

        star_ra = []
        star_dec = []
        star_mag = []

        for row in blend_input:
            block_center_x = row['xBlock']*1000 - 1500
            block_center_y = row['yBlock']*1000 - 1500
            if row['NumStars'] == 1:
                ra = (0.2*block_center_x) / 3600  # Convert to degrees
                dec = (0.2*block_center_y) / 3600
                star_ra.append(ra)
                star_dec.append(dec)
                star_mag.append(row['MagFaint'])
            else:
                star_x = []
                star_y = []
                row_mags = -1.*np.arange(row['NumStars'])*row['MagDiff']
                row_mags += row['MagFaint']
                rand_state.shuffle(row_mags)
                star_x.append(0.)
                star_y.append(0.)
                shift_angle = 0.
                for star_add in range(1, np.int(row['NumStars'])):
                    shift_angle = rand_state.uniform(-np.pi/2., np.pi/2.) + shift_angle
                    print(shift_angle, star_x[-1], star_y[-1])
                    star_x.append(row['PixSep']*np.cos(shift_angle)+star_x[-1])
                    star_y.append(row['PixSep']*np.sin(shift_angle)+star_y[-1])
                star_mean_x = np.mean(star_x)
                star_mean_y = np.mean(star_y)
                star_x = np.array(star_x) - star_mean_x + block_center_x
                star_y = np.array(star_y) - star_mean_y + block_center_y
                ra = (star_x*.2) / 3600
                dec = (star_y*.2) / 3600
                for i in range(np.int(row['NumStars'])):
                    star_ra.append(ra[i])
                    star_dec.append(dec[i])
                    star_mag.append(row_mags[i])

        star_ra = np.array(star_ra)
        star_dec = np.array(star_dec)
        star_ra += raInDegList[4]  # Use center CCD for now
        star_dec += declInDegList[4]
    
        for i in range(len(star_ra)):
            skySim.addStarByRaDecInDeg(i, star_ra[i],
                                       star_dec[i], star_mag[i])

        return skySim


if __name__ == "__main__":

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Create a Star Catalog for comcamCloseLoop simulation.")
    parser.add_argument("--inputFile", type=str,
                        help="input star layout file")
    parser.add_argument("--output", type=str, default="",
                        help="output catalog filename")
    args = parser.parse_args()

    if (args.output == ""):
        outputDir = getAoclcOutputPath()
        outputFilePath = os.path.join(outputDir, "starCat.txt")
    else:
        outputFilePath = args.output
        #outputFilePath = os.path.join(outputFilePath, "starCat.txt")

    createCat = createBlendTestCatalog()
    createCat.createPhosimCatalog(args.inputFile, outputFilePath)
