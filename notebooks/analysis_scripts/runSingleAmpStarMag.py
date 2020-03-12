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
    parser.add_argument("--opd", default=True, action='store_false')
    parser.add_argument("--defocalImg", default=True, action='store_false')
    parser.add_argument("--flats", default=True, action='store_false')
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()
    testLabel = 'mag'


    # settings for simulation
    numPro = 90 # number of processors setting in phosimCmptSetting.yaml 
    iterNum  = 1 # number of iterations 
    numFields = 9 # 9 for all CCDs, 3 is the minimum 

    # We create a PhoSim catalog with 2 stars with magVal brightness,
    # with varying separation in degrees     
    numStars  = 2 
    args.opd = False
    args.flats = False

    # since at such small separations the donuts overlap, we need to 
    # turn on the deblending ....
    doDeblending = True 

    # we want to save postage stamps..
    postageImg = True 

    # dir from /analysis_scripts/ level... 
    # - that's where we save the results
    topDir = 'results_after_centroid_update'
    expDir = 'singleAmpMag' # name of the experiment dir 

    # dir from /analysis_scripts/ level 
    # - that's where we copy the flats, opd from:
    copyDir = 'results_before_centroid_update/singleAmpSep'

    # the opd and wfs are stored here 
    os.environ["closeLoopTestDir"] = os.path.join(topDir, expDir) 

    starSep = 5 # % of amplifier span, i.e. 1.55 R_D - 
                #  then they're definitely overlappig
    magVal = 16 # brightness of primary 

    # iterate the brightness of secondary : 
    for mag in [16,15,14,13,12,11,10]:
        magList = [magVal, mag]
        print('\nStarting ccLoop for %d %% sep, using stars\
                with %d and %d mag'%(starSep, magVal, mag))

        outputDir = os.path.join(topDir,expDir,'mag_%d' % mag)
        print(outputDir)
        if (not os.path.exists(outputDir)):
            os.makedirs(outputDir)
       
        # re-use the flats and calibs : copy from the sep_10  ...
        if not args.opd and not args.flats : 
            print('Copying content of /sep_10/ to re-use the flats and OPD files...')
            argString = '-a '+os.path.join(copyDir,'sep_10')+'/. '+\
                         outputDir+'/'
            runProgram("cp", argstring=argString)

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
        skyFilePath = os.path.join(topDir,expDir,'starCat_%d.txt'%mag)
        print('The skyFilePath is %s'%skyFilePath)
        createCat = createPhosimCatalog()
        createCat.createPhosimCatalog(numStars, starSep, magList,skyFilePath, 
                  numFields=numFields)
        
        # initialize the baseComcamLoop.py Class 
        ccLoop = comcamLoop() 
        testName  = '%s.%d' % (testLabel, mag)
        print('The testName is %s'%testName)
        ccLoop.main(phosimDir, numPro, iterNum, outputDir, testName, 
                    isEimg=False,  genOpd=args.opd, genDefocalImg=args.defocalImg, 
                    genFlats=args.flats, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
                    doDeblending=doDeblending, postageImg=postageImg )
        print('Done running ccLoop for mag %d\n\n'%mag)
 
        # # Once the necessary data is created we don't need to recreate on every iteration
        args.opd = False
        args.flats = False
        #args.defocalImg = False
