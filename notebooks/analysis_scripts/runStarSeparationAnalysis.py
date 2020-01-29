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
    parser.add_argument("--testLabel", type=str, default="sep")
    parser.add_argument("--testOutput", type=str, default="")
    parser.add_argument("--skyFile", type=str, default="starCat.txt")
    parser.add_argument("--raShift", type=float, default=0.0)
    parser.add_argument("--decShift", type=float, default=0.0)
    parser.add_argument("--magVal", type=float, default=15)
    # parser.add_argument("--opd", default=True, action='store_false')
    # parser.add_argument("--defocalImg", default=True, action='store_false')
    # parser.add_argument("--flats", default=True, action='store_false')
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()
    outputDir = getAoclcOutputPath()
    testLabel = args.testLabel
    skyFilePath = args.skyFile


   

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir
 

    # We create a PhoSim catalog with 2 stars with magVal brightness,
    # with varying separation in degrees     
    for starSep in np.arange(0.01, 0.2347, 0.025 ) : 

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
        

        # read the star brightness
        magVal   = args.magVal
        
        createCat.createPhosimCatalog(2, starSep, [magVal,magVal], raShift, decShift,
                                      skyFilePath)

        ccLoop = comcamLoop() 
        numPro = 8 # number of processors setting in phosimCmptSetting.yaml 
        iterNum  = 1 # number of iterations 
        ccLoop.main(phosimDir, numPro, 1, outputDir, '%s.%.1f' % (testLabel, starSep), 
                    isEimg=False, genOpd=True, genDefocalImg=True, 
                    genFlats=True, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05)

        # # Once the necessary data is created we don't need to recreate on every iteration
        args.opd = False
        args.flats = False
        args.defocalImg = False
