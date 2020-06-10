import os
import argparse
import numpy as np
from baseComcamLoop import baseComcamLoop as comcamLoop
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram

# Load directory paths
phosimDir = getPhoSimPath()
testLabel = 'sep'

# set the star brightness
magVal   = 16

# settings for simulation
numPro = 20 # number of processors setting in phosimCmptSetting.yaml 
iterNum  = 1 # number of iterations 
numFields = 9 # 9 for all CCDs,  3 to get the result quicker... 

# We create a PhoSim catalog with 2 stars with magVal brightness,
# with varying separation 
numStars  = 2 
opd = True # initially we calculate opd and flats .. 
flats = True
defocalImg = True

# since at such small separations the donuts overlap, we need to 
# turn on the deblending ....
doDeblending = True 

# do we want to save postage stamps in WepController.py  getDonutMap() ..
postageImg = True 

# dir from /analysis_scripts/ level... 
topDir = 'results_after_centroid_update'
expDir = 'singleAmpSepNew' # name of the experiment dir 

# the opd and wfs are stored here 
os.environ["closeLoopTestDir"] = os.path.join(topDir, expDir) 

# Run the same OPD simulation for varying star separation ... 

# 
magList = np.ones(numStars) * magVal

# separation defined as % of  ra amplifier span (eg. 10%, 50%...)
# we start from largest separation and move closer 
for starSep in np.arange(1,11)[::-1] : 
    print('\nStarting ccLoop for separation equal to %d percent of \
        amplifier ra span '%starSep)

    outputDir = os.path.join(topDir,expDir,'sep_%d' % starSep)
    print(outputDir)
    if (not os.path.exists(outputDir)):
        os.makedirs(outputDir)

    # first time running the loop, opd=True, flats=True,
    # then it is set to opd=False, flats=False, to
    # re-use the flats and calibs : copy from the sep_10  ...

    if not opd and not flats : 
        print('Copying content of /sep_10/ to re-use the flats and OPD files...')
        argString = '-a '+os.path.join(topDir,expDir,'sep_10')+'/. '+\
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

    # Clobber everything
    if opd is True:
        print('We will make new OPD files in this run')
        _eraseFolderContent(outputDir)

    else: 

        if flats is True: # Clobber just calibs  
            print('We will make new flats in this run')
            _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
            _eraseFolderContent(os.path.join(outputDir, 'input'))  

        if defocalImg is True: # Clobber only defocal images 
            print('We will make new defocal images in this run ')
            intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
            extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
            if os.path.exists(intraPath):
                _eraseFolderContent(intraPath)
            if os.path.exists(extraPath):
                _eraseFolderContent(extraPath)


    # make the star catalog 
    skyFilePath = os.path.join(topDir,expDir,'starCat_%d.txt'%starSep)
    print(skyFilePath)
    createCat = createPhosimCatalog()
    createCat.createPhosimCatalog(numStars, starSep, magList,skyFilePath,numFields=numFields,
                                 addStarsInAmp=True)

    # #initialize the baseComcamLoop.py Class 
    ccLoop = comcamLoop() 

    print('For separation of %d %%,  the outputDir is %s'%(starSep,outputDir))
    testName = '%s.%d' % (testLabel, starSep)
    ccLoop.main(phosimDir, numPro, iterNum, outputDir, testName, 
                isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, 
                genFlats=flats, useMinDofIdx=False,
                inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
                doDeblending=doDeblending, postageImg=postageImg )
    print('Done running ccLoop for this separation\n\n')

    # # Once the necessary data is created we don't need to recreate on every iteration
    opd = False
    flats = False