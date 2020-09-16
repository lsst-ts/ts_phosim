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
numPro = 15 # number of processors setting in phosimCmptSetting.yaml 
iterNum  = 1 # number of iterations 

# we evaluate OPD for LsstFamCam
opd = True 

# we simulate the calibration products for LsstCam to do the ISR 
flats = True

# we simulate defocal corner LsstCam images 
defocalImg = True 

# no need to deblend single stars 
doDeblending = True 

# we want to save postage stamps..
postageImg = True 

# dir from /analysis_scripts/ level... 
# - that's where we save the results
topDir = 'results_wfs'
expDir = 'arrowStars_2020_24' # name of the experiment dir 

outputDir = os.path.join(topDir,expDir)

# Use Te-Wei's catalog
skyFilePath = os.path.join(topDir, 'skyWfsArrow.txt')

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
   

# initialize the baseWfsWep.py Class 
wfsLoop = baseWfsWep() 
testName  = '%s' % (testLabel)
print('The testName is %s'%testName)

wfsLoop.main(phosimDir, numPro, iterNum, outputDir, testName, 
            isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, 
            genFlats=flats, useMinDofIdx=False,
            inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
            doDeblending=doDeblending, postageImg=postageImg,
            phosimRepackagerKeepOriginal = True )

print('Done running wfsLoop \n\n')


