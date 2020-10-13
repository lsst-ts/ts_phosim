# old version of the star separation test code 

# import sys
# path_to_comcamCloseLoop = '/astro/store/epyc/projects/lsst_comm/ts_phosim/bin.src/'
# sys.path.append(path_to_comcamCloseLoop)
import os
import argparse
import numpy as np
from baseComcamLoop import baseComcamLoop as comcamLoop
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram

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
    numPro = 24 # number of processors setting in phosimCmptSetting.yaml 
    iterNum  = 1 # number of iterations 
    numFields = 3 # 9 for all CCDs,  3 to get the result quicker... 

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir
 

    # We create a PhoSim catalog with 2 stars with magVal brightness,
    # with varying separation in degrees     
    numStars  = 2 
    args.opd = False
    args.flats = False

    # since at such small separations the donuts overlap, we need to 
    # turn on the deblending ....
    doDeblending = True 

    magList = np.ones(numStars) * magVal
    for starSep in [1] : 
        print('\nStarting ccLoop for separation equal to %d percent of amplifier ra span '%starSep)

        outputDir = 'singleAmp/sep_%d' % starSep
        print(outputDir)
        if (not os.path.exists(outputDir)):
            os.makedirs(outputDir)
        
        # re-use the flats and calibs : copy from the sep_10  ...

        if not args.opd and not args.flats : 
            #print('Copying content of /sep_10/ to re-use the flats and OPD files...')
            #argString = '-a singleAmp/sep_10/. '+ outputDir+'/'
            #runProgram("cp", argstring=argString)

            # ensure that input/raw and input/rerun are empty 
            print('Deleting content of input/raw/ and input/rerun/')
            _eraseFolderContent(os.path.join(outputDir, 'input','raw'))
            _eraseFolderContent(os.path.join(outputDir, 'input','rerun'))

            # remove files that are remade
            argString = os.path.join(outputDir, 'input')
            runProgram("rm", argstring=argString+'/isr*')
            runProgram("rm", argstring=argString+'/registry*')
            runProgram("rm", argstring=argString+'/_mappe*')

        # Clobber
        if args.opd is True:
            print('We will make new OPD files in this run')
            _eraseFolderContent(outputDir)
        else:
            if args.flats is True:
                print('We will make new flats in this run')
                _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
                _eraseFolderContent(os.path.join(outputDir, 'input'))     
            if args.defocalImg is True:
                print('We will make new defocal images in this run ')
                intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
                extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
                if os.path.exists(intraPath):
                    _eraseFolderContent(intraPath)
                if os.path.exists(extraPath):
                    _eraseFolderContent(extraPath)
        
       
        # make the star catalog 
        skyFilePath = 'singleAmp/starCat_%d.txt'%starSep
        print(skyFilePath)
        createCat = createPhosimCatalog()
        createCat.createPhosimCatalog(numStars, starSep, magList,skyFilePath,numFields=numFields)
        
        # #initialize the baseComcamLoop.py Class 
        ccLoop = comcamLoop() 
        
        print('For separation of %d %%,  the outputDir is %s'%(starSep,outputDir))
        ccLoop.main(phosimDir, numPro, iterNum, outputDir, '%s.%d' % (testLabel, starSep), 
                    isEimg=False,  genOpd=args.opd, genDefocalImg=args.defocalImg, 
                    genFlats=args.flats, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
                    doDeblending=doDeblending)
        print('Done running ccLoop for this separation\n\n')
 
        # # Once the necessary data is created we don't need to recreate on every iteration
        #args.opd = False
        #args.flats = False
        #args.defocalImg = False
