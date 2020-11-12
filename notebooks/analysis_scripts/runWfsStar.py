import os
import numpy as np
from baseLsstCamLoop import baseLsstCamLoop
from baseComcamLoop import _eraseFolderContent
from baseAOSLoop import baseAOSLoop
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram

# Load directory paths
phosimDir = getPhoSimPath()
testLabel = 'wfs'

# settings for simulation
numPro = 20 # number of processors setting in phosimCmptSetting.yaml 
iterNum  = 5 # number of iterations 

# we evaluate OPD for 31 locations of LsstFamCam
opd = True 

# we simulate the calibration products for LsstCam to do the ISR 
flats = True

# we simulate defocal corner LsstCam images 
defocalImg = True 

# no need to deblend single stars 
doDeblending = False 

# we want to save postage stamps..
postageImg = True 

# dir from /analysis_scripts/ level... 
# - that's where we save the results
topDir = 'results_wfs'
expDir = 'arrowStarsLetters_2020_44' # name of the experiment dir 

outputDir = os.path.join(topDir,expDir)

# the  catalog to use 
skyFilePath = os.path.join(topDir, 'skyWfsArrowLetter.txt')

# the opd and wfs are stored here 
os.environ["closeLoopTestDir"] = outputDir


print('\nStarting wfsWep, the outputDir is %s '%outputDir)

if (not os.path.exists(outputDir)):
    os.makedirs(outputDir)

if flats is True:
    print('We will make new flats in this run')
    flatDir  = os.path.join(outputDir, 'fake_flats')
    inputDir = os.path.join(outputDir, 'input')
    if os.path.exists(flatDir):  _eraseFolderContent(flatDir)
    if os.path.exists(inputDir): _eraseFolderContent(inputDir)    
else:
    ingestedFlats = os.path.join(outputDir,'input','flat') 
    if os.path.exists(ingestedFlats):
        print('Assuming there are already ingested flats in input/flat/')
    else:
        raise RuntimeError("No ingested flats existing ...")    

if defocalImg is True:
    print('We will make new defocal images in this run ')
    intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
    extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
    if os.path.exists(intraPath):  _eraseFolderContent(intraPath)
    if os.path.exists(extraPath):  _eraseFolderContent(extraPath)   
   
opdCmd  = 'opdDefault.cmd'
starCmd = 'starDefault.cmd'

# initialize the baseWfsWep.py Class 
wfsLoop = baseAOSLoop() 

wfsLoop.main(phosimDir, numPro, iterNum, outputDir, testLabel, 
            isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, 
            genFlats=flats, useMinDofIdx=False,
            inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
            doDeblending=doDeblending, postageImg=postageImg,
            phosimRepackagerKeepOriginal = True, opdCmdSettingsFile=opdCmd, 
            starCmdSettingsFile=starCmd, selectSensors='lsstcam',
            expWcs = False, noPerturbations=False )

print('\nDone running lsstCamLoop \n\n')


