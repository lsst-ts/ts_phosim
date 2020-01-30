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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--testLabel", type=str, default="mag")
    parser.add_argument("--testOutput", type=str, default="")
    parser.add_argument("--skyFile", type=str, default="starCat.txt")
    parser.add_argument("--raShift", type=float, default=0.0)
    parser.add_argument("--decShift", type=float, default=0.0)
    # parser.add_argument("--opd", default=True, action='store_false')
    # parser.add_argument("--defocalImg", default=True, action='store_false')
    # parser.add_argument("--flats", default=True, action='store_false')
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()
    outputDir = getAoclcOutputPath()
    testLabel = args.testLabel
    skyFilePath = args.skyFile
    genFlats = True

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir

    for magVal in np.arange(10.0, 16.1, 0.5):

        # # Clobber
        # if args.opd is True:
        _eraseFolderContent(outputDir)
        # else:
        #     if args.flats is True:
        #         _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
        #         _eraseFolderContent(os.path.join(outputDir, 'input'))     
        #     if args.defocalImg is True:
        #         _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'intra'))
        #         _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'extra'))

        createCat = createPhosimCatalog()
        raShift = (args.raShift * .2) / 3600 # Convert to degrees
        decShift = (args.decShift * .2) / 3600 # Convert to degrees
        createCat.createPhosimCatalog(1, 0, [magVal], raShift, decShift,
                                      skyFilePath, numFields=3)

        ccLoop = comcamLoop()
        ccLoop.main(phosimDir, 8, 1, outputDir, '%s.%.1f' % (testLabel, magVal), 
                    isEimg=False, genOpd=True, genDefocalImg=True, 
                    genFlats=genFlats, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05)

        # # Once the necessary data is created we don't need to recreate on every iteration
        # args.opd = False
        # genFlats = False
        # args.defocalImg = False
