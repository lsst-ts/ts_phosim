import os
import numpy as np
from baseComcamLoop import baseComcamLoop as comcamLoop
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram

# initial setting whether  to calculate opd, etc. 
flats = False
opd = False
defocalImg = True
justWfs = False # switch to only re-do wfs,  

simulatedSensors  = 'comcam' # 'wfs'

#not making or copying flats, opd, defocalImg...

# just to be consistent...
if justWfs:
    flats=  False
    opd = False
    defocalImg = False # don't re-generate opd
    # if just calculating wfs ... 

# Load directory paths
phosimDir = getPhoSimPath()
testLabel = 'gaia'

# settings for simulation
numPro = 20 # of 90 : the number of processors to use  in phosimCmptSetting.yaml 
iterNum  = 1 # number of iterations 
numFields = 9 # 9 for all CCDs,  3 to get the result quicker... 


# We use as PhoSim catalog a GAIA star catalog ...  
 
# since at such small separations the donuts overlap, we need to 
# turn on the deblending ....
doDeblending = True 


# change the camDimOffset  setting in ts_wep/policy
# to not select stars that are too close to the edge ...
camDimOffset = -50


# we want to save postage stamps - the 
# default directory is 
# outputDir + '/postage'
postageImg = True 

# dir from /analysis_scripts/ level
# - that's where we save the results:
topDir = 'results_gaia'
expDir = 'gMagGt11_R22_' # name of the experiment dir 

# dir from /analysis_scripts/ level 
# - that's where we copy the flats, opd from:
#copyDir = 'results_before_centroid_update/singleAmpSep/sep_10'
copyDir = 'results_gaia/gMagGt11'

# the opd and wfs are stored here 
os.environ["closeLoopTestDir"] = os.path.join(topDir, expDir) 

print('\nStarting ccLoop for GAIA catalog')

outputDir = os.path.join(topDir,expDir)
print('The outputDir is %s'%outputDir)
if (not os.path.exists(outputDir)):
    os.makedirs(outputDir)



if justWfs:
    print('Skipping  OPD and flat copying, not making new defocalImg')
    print('Just calculating WFS ')
 # only do all that if not trying to just rerun the WFS ... 
if not justWfs :  
    # re-use the flats and calibs,
    # since these have nothing to do with the actual stars ...
    # -  copy from the sep_10  ...
    if not opd and not flats : 
        print('Copying content of %s to re-use the flats and OPD files...'%copyDir)
        
        # first copy the old results...
        argString = '-a '+copyDir+'/. '+outputDir+'/'
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
    if opd is True:
        print('We will make new OPD files in this run')
        _eraseFolderContent(outputDir)
    else:
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


# read the star catalog from file ... 
# it conforms to the format expected by PhoSim 
# mv /data/epyc/users/suberlak/starCatGAIA.txt analysis_scripts/results_gaia
if simulatedSensors is 'comcam':
    skyFilePath = os.path.join(topDir,'starCatGAIA_gt11.txt')
    onlyComcam  = True
    onlyWfsSensors = False

if simulatedSensors is 'wfs':
    skyFilePath  = os.path.join(topDir, 'starCatGAIA_gt13_wfs.txt')
    onlyComcam = False
    onlyWfsSensors = True 

print('Using sky catalog from %s'%skyFilePath)
       
# set the opd.cmd and star.cmd files ...
opdCmd  = 'opdQuickBackground.cmd'
comcamCmd = 'starQuickBackground.cmd'

print('For PhoSim using /policy/cmdFile/%s and %s'%(opdCmd,comcamCmd))

# initialize the baseComcamLoop.py Class 
ccLoop = comcamLoop() 
ccLoop.main(phosimDir, numPro, iterNum, outputDir, testLabel, 
            isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, 
            genFlats=flats, useMinDofIdx=False,
            inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
            doDeblending=doDeblending, camDimOffset = camDimOffset, 
            postageImg=postageImg ,
            opdCmdSettingsFile=opdCmd,
            comcamCmdSettingsFile=comcamCmd,onlyComcam = onlyComcam, 
            onlyWfsSensors  = onlyWfsSensors)
print('Done running ccLoop for GAIA \n\n')
