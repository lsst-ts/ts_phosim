# the version used for AOS_singleAmpMag.ipynb
import os
import argparse
import numpy as np
from baseComcamLoop import baseComcamLoop as comcamLoop
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram

# parse args so that the program can be run from within the notebook,
# with arguments stored in the notebook rather than directly in this file 

parser = argparse.ArgumentParser(
    description='Run ComCam AOS closed loop for varying star magnitude at fixed separation')
parser.add_argument("--testLabel", type=str, default="mag",
                    help='test label')
parser.add_argument("--magPrimary", type=float, default=16, 
                    help='magnitude for the primary (fixed)')
parser.add_argument("--starSep", type=float, default=5,
                   help='star separation, in percentage of ra span of the amplifier - \
                   for the default 5 %%, the donuts overlap.')
parser.add_argument("--numPro", type=int, default=10, 
                    help='number of processors setting in phosimCmptSetting.yaml')
parser.add_argument("--iterNum", type=int, default=1, 
                    help='number of iterations ')
parser.add_argument("--numFields", type=int, default=9, 
                    help='number of CCDs, 3 is the minimum for ts_wep')
parser.add_argument("--numStars", type=int, default=2, 
                    help='number of stars to simulate')
parser.add_argument("--opd", default=True, action='store_true',
                   help='whether to calculate the truth==the optical path difference. \
                   If False, it seeks to copy OPD from copyDir')
parser.add_argument("--flats", default=True, action='store_true',
                   help='whether to calculate calibration files. \
                   If False, it seeks to copy calibs from copyDir')
parser.add_argument("--defocalImg", default=True, action='store_true',
                   help='whether to simulate new defocal images')
parser.add_argument("--doDeblending", default=True, action='store_true',
                   help='whether to use deblending code in ts_wep')
parser.add_argument("--postageImg", default=True, action='store_true',
                   help='whether to store postage stamp images \
                   in /postage/, prepared in WepController.py')
parser.add_argument("--topDir", type=str, default='results_after_centroid_update',
                   help='the access directory from the level of running this code')
parser.add_argument("--expDir", type=str, default='singleAmpMagNew',
                   help='the experiment directory within topDir, \
                   the closeLoopTestDir=os.path.join(topDir,expDir)')
parser.add_argument("--copyDir", type=str, default='mag_16',
                   help='the directory within topDir, \
                   from which we copy calibs and OPD if flats=False and opd=False.\
                   By default, the widest separation which is run first in the \
                   loop over star separations.')

args = parser.parse_args()

# Load directory paths
phosimDir = getPhoSimPath()
testLabel = args.testLabel

# settings for simulation
numPro = args.numPro # number of processors setting in phosimCmptSetting.yaml 
iterNum  = args.iterNum # number of iterations 
numFields = args.numFields # 9 for all CCDs,  3 to get the result quicker... 

# We create a PhoSim catalog with 2 stars with fixed separation,
# and brightness of the primary, 
# varying the brighness of the secondary brightness, 
numStars  = args.numStars  
opd = args.opd # initially we calculate opd and flats .. 
flats = args.flats
defocalImg = args.defocalImg

# since at such small separations the donuts overlap, we need to 
# turn on the deblending ....
doDeblending = args.doDeblending 

# do we want to save postage stamps in WepController.py  getDonutMap() ..
postageImg = args.postageImg 

topDir = args.topDir # dir from /analysis_scripts/ level... 
expDir = args.expDir # name of the experiment dir 
copyDir = args.copyDir # name of dir to copy calibs and OPD from 

# the opd and wfs are stored here 
os.environ["closeLoopTestDir"] = os.path.join(topDir, expDir) 

starSep = args.starSep  # % of amplifier span, i.e. 1.55 R_D - 
                #  then they're definitely overlappig
magPrimary = args.magPrimary # brightness of primary 

# iterate the brightness of secondary : 
for magSecondary in [16,15,14,13,12,11,10]:
    magList = [magPrimary, magSecondary]
    print('\nStarting ccLoop for %d %% sep, using stars\
            with %d and %d mag'%(starSep, magPrimary, magSecondary))

    outputDir = os.path.join(topDir,expDir,'mag_%d' % magSecondary)
    print(outputDir)
    if (not os.path.exists(outputDir)):
        os.makedirs(outputDir)
       
    
    # first time running the loop, opd=True, flats=True,
    # then it is set to opd=False, flats=False, to
    # re-use the flats and calibs : copy from the copyDir

    if not opd and not flats : 
        print('Copying content of %s to re-use the flats and OPD files...'%copyDir)

        argString = '-a '+os.path.join(topDir,expDir,copyDir)+'/. '+\
                     outputDir+'/'
        runProgram("cp", argstring=argString)

        # ensure that input/raw and input/rerun are empty 
        print('Deleting content of input/raw/ and input/rerun/')
        _eraseFolderContent(os.path.join(outputDir, 'input','raw'))
        _eraseFolderContent(os.path.join(outputDir, 'input','rerun'))

        # deleting postage images 
        print('Deleting content of /postage/ ')
        _eraseFolderContent(os.path.join(outputDir,'postage'))

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
    skyFilePath = os.path.join(topDir,expDir,'starCat_%d.txt'%magSecondary)
    print('The skyFilePath is %s'%skyFilePath)
    createCat = createPhosimCatalog()
    createCat.createPhosimCatalog(outputFilePath=skyFilePath,
                              numStars=numStars, numFields=numFields,
                              addStarsInAmp=True, starSep=starSep, 
                              magList=magList, selectSensors='comcam')

    # initialize the baseComcamLoop.py Class 
    ccLoop = comcamLoop() 
    testName  = '%s.%d' % (testLabel, magSecondary)
    print('For magPrimary %.2f and magSecondary %.2f testName is %s'%(magPrimary, magSecondary, testName))
    ccLoop.main(phosimDir=phosimDir, numPro=numPro, iterNum=iterNum, 
                baseOutputDir=outputDir, testName=testName,
                isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, 
                genFlats=flats, useMinDofIdx=False,
                inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
                doDeblending=doDeblending, postageImg=postageImg )
    print('Done running ccLoop for this magnitude of secondary: %.2f\n\n'%magSecondary)

     # # Once the necessary data is created we don't need to recreate on every iteration
    opd = False
    flats = False
