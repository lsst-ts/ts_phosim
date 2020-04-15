# import sys
# path_to_comcamCloseLoop = '/astro/store/epyc/projects/lsst_comm/ts_phosim/bin.src/'
# sys.path.append(path_to_comcamCloseLoop)
import os
import numpy as np
from baseWfsWep import baseWfsWep
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram

# Load directory paths
phosimDir = getPhoSimPath()
testLabel = 'wfs'

# settings for simulation
numPro = 10 # number of processors setting in phosimCmptSetting.yaml 
iterNum  = 1 # number of iterations 
numFields = 9 # 9 for all CCDs, 3 is the minimum 

# We create a PhoSim catalog with 1 star on each WFS sensor
# we do not attempt to calculate OPD...
numStars  = 2 
opd = False 

# But we can calculate the calibration products to do the ISR 
flats = False


# And we definitely want to calculate defocal image to run WEPCalc ....
defocalImg = True 

# no need to deblend single stars 
doDeblending = False 

# we want to save postage stamps..
postageImg = True 

# dir from /analysis_scripts/ level... 
# - that's where we save the results
topDir = 'results_wfs'
expDir = 'singleStar' # name of the experiment dir 

# dir from /analysis_scripts/ level 
# - that's where we copy the flats, opd from:
#copyDir = 'results_before_centroid_update/singleAmpSep'

outputDir = os.path.join(topDir,expDir)


# the opd and wfs are stored here 
os.environ["closeLoopTestDir"] = outputDir

starMag = 16 # brightness of primary 

print('\nStarting wfsWep, the outputDir is %s '%outputDir)

if (not os.path.exists(outputDir)):
    os.makedirs(outputDir)

# re-use the flats and calibs : copy from the sep_10  ...
# if not args.opd and not args.flats : 
#     print('Copying content of /sep_10/ to re-use the flats and OPD files...')
#     argString = '-a '+os.path.join(copyDir,'sep_10')+'/. '+\
#                  outputDir+'/'
#     runProgram("cp", argstring=argString)

#     # ensure that input/raw and input/rerun are empty 
#     print('Deleting content of input/raw/ and input/rerun/')
#     _eraseFolderContent(os.path.join(outputDir, 'input','raw'))
#     _eraseFolderContent(os.path.join(outputDir, 'input','rerun'))

#     # remove files that are remade
#     argString = os.path.join(outputDir, 'input')
#     runProgram("rm", argstring=argString+'/isr*')
#     runProgram("rm", argstring=argString+'/registry*')
#     runProgram("rm", argstring=argString+'/_mappe*')

# Clobber
# if  opd is True:
#     print('We will make new OPD files in this run')
#     _eraseFolderContent(outputDir)
# else:
if flats is True:
    print('We will make new flats in this run')
    _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
    _eraseFolderContent(os.path.join(outputDir, 'input'))     
if defocalImg is True:
    print('We will make new defocal images in this run ')
    intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
    extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
    if os.path.exists(intraPath):
        _eraseFolderContent(intraPath)
    if os.path.exists(extraPath):
        _eraseFolderContent(extraPath)   
   
# make the star catalog 
skyFilePath = os.path.join(outputDir, 'starCat.txt')
print('The skyFilePath is %s'%skyFilePath)
createCat = createPhosimCatalog()
createCat.createPhosimCatalog(outputFilePath=skyFilePath, starMag=starMag, 
                              selectSensors='wfs')


# initialize the baseWfsWep.py Class 
wfsLoop = baseWfsWep() 
testName  = '%s.%d' % (testLabel, starMag)
print('The testName is %s'%testName)

wfsLoop.main(phosimDir, numPro, iterNum, outputDir, testName, 
            isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, 
            genFlats=flats, useMinDofIdx=False,
            inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
            doDeblending=doDeblending, postageImg=postageImg )
print('Done running wfsLoop for mag %d\n\n'%starMag)


