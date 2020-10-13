# OLD version of star separation analusis (first try ... )

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
    #testOutputDir = '/data/epyc/users/suberlak/Commissioning/aos/aoclc_analysis/'
    parser.add_argument("--testLabel", type=str, default="sep")
    parser.add_argument("--testOutput", type=str, default="")
    parser.add_argument("--skyFile", type=str, default="starCat.txt")
    parser.add_argument("--raShift", type=float, default=150.0)
    parser.add_argument("--decShift", type=float, default=150.0)
    parser.add_argument("--magVal", type=float, default=16)
    parser.add_argument("--opd", default=True, action='store_false')
    parser.add_argument("--defocalImg", default=True, action='store_false')
    parser.add_argument("--flats", default=True, action='store_false')
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()
    #outputDir = getAoclcOutputPath()
    testLabel = args.testLabel
    #skyFilePath = args.skyFile
   

    # read the pixel shift off center 
    raShift = (args.raShift * .02) / 3600 # Convert to degrees
    decShift = (args.decShift * .02) / 3600 # Convert to degrees

    # read the star brightness
    magVal   = args.magVal


    # settings for simulation
    numPro = 60 # number of processors setting in phosimCmptSetting.yaml 
    iterNum  = 1 # number of iterations 
    numFields = 9 # 9 for all CCDs,  3 to get the result quicker... 

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir
 

    # We create a PhoSim catalog with 2 stars with magVal brightness,
    # with varying separation in degrees     
    starNum  = 2 
    for starSep in np.arange(0.01, 0.02, 0.025 ): 
        print('\nStarting ccLoop for separation %.3f'%starSep)

        outputDir = 'output/sep_%.3f' % starSep
        print(outputDir)
        if (not os.path.exists(outputDir)):
            os.makedirs(outputDir)


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
        
       
        # make the star caatalog 
        skyFilePath = 'output/starCat_%.3f.txt'%starSep
        print(skyFilePath)
        createCat = createPhosimCatalog()
        createCat.createPhosimCatalog(starNum, starSep, [magVal,magVal], raShift, decShift,
                                      skyFilePath,numFields=numFields)
        
        # #initialize the baseComcamLoop.py Class 
        ccLoop = comcamLoop() 
        
        print('For starSep%f, the outputDir is %s'%(starSep,outputDir))
        ccLoop.main(phosimDir, numPro, iterNum, outputDir, '%s.%.3f' % (testLabel, starSep), 
                    isEimg=False,  genOpd=args.opd, genDefocalImg=args.defocalImg, 
                    genFlats=args.flats, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05)
        print('Done running ccLoop for this separation\n\n')
 
        # # Once the necessary data is created we don't need to recreate on every iteration
        args.opd = False
        args.flats = False
        args.defocalImg = False
