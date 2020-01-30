# import sys
# path_to_comcamCloseLoop = '/astro/store/epyc/projects/lsst_comm/ts_phosim/bin.src/'
# sys.path.append(path_to_comcamCloseLoop)
import os
import argparse
import numpy as np
from baseComcamLoop import baseComcamLoop as comcamLoop
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalog import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir

class centroidLoop(comcamLoop):

    def __init__(self, centroidInfoFileName):

        self.centroidInfoFileName = centroidInfoFileName

    def _outputSkyInfo(self, outputDir, skyInfoFileName,
                       skySim, wepCalc):

        outputSkyInfoFilePath = os.path.join(outputDir, skyInfoFileName)
        skySim.exportSkyToFile(outputSkyInfoFilePath)
        wepCalc.setSkyFile(self.centroidInfoFileName)

        return skySim, wepCalc


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--testLabel", type=str, default="centroid")
    parser.add_argument("--testOutput", type=str, default="")
    parser.add_argument("--skyFile", type=str, default="starCat.txt")
    parser.add_argument("--centroidFile", type=str, default="perturbStarCat.txt")
    parser.add_argument("--raShift", type=float, default=0.0)
    parser.add_argument("--decShift", type=float, default=0.0)
    parser.add_argument("--opd", default=True, action='store_false')
    parser.add_argument("--defocalImg", default=True, action='store_false')
    parser.add_argument("--flats", default=True, action='store_false')
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()
    outputDir = getAoclcOutputPath()
    testLabel = args.testLabel
    skyFilePath = args.skyFile
    centroidSkyFilePath = args.centroidFile

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir

    for offset in range(-40, 41, 10):

        # Clobber
        if args.opd is True:
            _eraseFolderContent(outputDir)
        else:
            if args.flats is True:
                _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
                _eraseFolderContent(os.path.join(outputDir, 'input'))     
            if args.defocalImg is True:
                _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'intra'))
                _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'extra'))

        createCat = createPhosimCatalog()
        raShift = (args.raShift * .2) / 3600 # Convert to degrees
        decShift = (args.decShift * .2) / 3600 # Convert to degrees
        createCat.createPhosimCatalog(1, 0, [15], raShift, decShift,
                                      skyFilePath)

        rand_state = np.random.RandomState(seed=(100+offset))
        offset_ang = 0 # rand_state.choice(np.arange(0, 360)) # 0 is y-direction (ra), 90 is x-direction (dec)
        pixelOffset = np.array([np.cos(np.radians(offset_ang))*offset,
                                np.sin(np.radians(offset_ang))*offset]) # ra is y direction as set up
        centroidOffset = (pixelOffset * .2) / 3600 # Convert to degrees
        createCat = createPhosimCatalog()
        createCat.createPhosimCatalog(1, 0, [15], raShift + centroidOffset[0], 
                                      decShift + centroidOffset[1],
                                      str(centroidSkyFilePath[:-3] + str(offset) + 'ang' + str(offset_ang) + '.txt'))

        ccLoop = centroidLoop(str(centroidSkyFilePath[:-3]+str(offset) + 'ang' + str(offset_ang) + '.txt'))
        ccLoop.main(phosimDir, 8, 1, outputDir, str(testLabel + '.' + str(offset) + 'ang'+ str(offset_ang)), 
                    isEimg=False, genOpd=args.opd, genDefocalImg=args.defocalImg, 
                    genFlats=args.flats, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05)

        # Once the necessary data is created we don't need to recreate on every iteration
        args.opd = False
        args.flats = False
        args.defocalImg = False
