import os
import datetime
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
focalImg  = True

copy = True


justWfs = False # switch to only re-do wfs,  
splitWfsByMag = False # whether to calculate the wfs for subsets of stars 
                     # based on magnitude ranges, or not ...
                     # the ranges are 11-16 mag, as the normal
                     # Filter.py g-mag range is 9.74-16.17,
                     # so even though there are stars in gaia DR2 catalog
                     # for the region centered on  boresight (ra,dec)=(0,0) 
                     # from 10 to 22 mag,  
                     # those fainter  than 16.17 are not used for WFS 
                     # calculation anyway ... 


selectSensors  = 'comcam' #  or 'wfs', or 'None' - if none, then 
# PhoSim will simulate all sensors that get any brightness on them, 
# even if it's just one photon - thus if there are any bright stars that 
# may be saturated, it will simulate these sensors, even if they 
# are beyond Raft22 , which is a problem  for the phosim_repackager,
# which is not expecting it ... 

#not making or copying flats, opd, defocalImg...

# just to be consistent...
if justWfs:
    print('We are in WFS-calculation only mode, re-using flats, opd, and defocal images ')
    flats=  False
    opd = False
    defocalImg = False # don't re-generate opd
    focalImg = False
    # if just calculating wfs ... 

# Load directory paths
phosimDir = getPhoSimPath()
testLabel = 'gaia'

# settings for simulation
numPro = 10 # of 90 : the number of processors to use  in phosimCmptSetting.yaml 
            # need to make sure this is many times less than the actual 
            # numProc b/c each process is massively parallel 
iterNum  = 1 # number of iterations 
numFields = 9 # 9 for all CCDs,  3 to get the result quicker... 


# We use as PhoSim catalog a GAIA star catalog ...  
 
# since at such small separations the donuts overlap, we need to 
# turn on the deblending ....
doDeblending = True 


# change the camDimOffset  setting in ts_wep/policy
# to not select stars that are too close to the edge ...
camDimOffset = -150


# we want to save postage stamps - the 
# default directory is 
# outputDir + '/postage'
postageImg = True 

# dir from /analysis_scripts/ level
# - that's where we save the results:
topDir = 'results_gaia'

# simulating four galactic locations:  
field = 'Baade' # 'med' lBaade'
catalogType = 'gt11' # 'full'

expDir = 'dr2_%s_%s'%(field,catalogType) # name of the experiment dir 
# 'dr2_med_gt11'
# 'dr2_low_gt11'
# 'dr2_baade_gt11'

# dir from /analysis_scripts/ level 
# - that's where we copy the flats, opd from:
#copyDir = 'results_before_centroid_update/singleAmpSep/sep_10'
copyDir = 'results_gaia/dr2_high_gt11'

print('\nStarting ccLoop for GAIA catalog')

outputDir = os.path.join(topDir,expDir)
print('The outputDir is %s'%outputDir)

# the opd and wfs are stored here 
os.environ["closeLoopTestDir"] = outputDir


if (not os.path.exists(outputDir)):
    os.makedirs(outputDir)


if justWfs:
    #print('Skipping  OPD and flat copying, not making new defocalImg')
    print('Just calculating WFS ')

if copy : 
    print('Copying entire content of %s ...'%copyDir)
    print('This includes calibs, OPD, defocal images, ingest, etc.')
    #print('NB: if phosim /LSST stack were recently updated,\
    #     you need to make sure that OPD, calibs were made and \
    #     ingested with the same version ! Otherwise you may see errors \
    #     with eg. phosim_repackager.py ')

    #first copy the old results...
    argString = '-a '+copyDir+'/. '+outputDir+'/'
    runProgram("cp", argstring=argString)
    
    # remove the ISR directory,
    # especially the input/rerun/run1/repositoryCfg.yaml , which 
    # contains the explicit absolute path 
    # eg. 
    # _root: /data/epyc/users/suberlak/Commissioning/aos/ts_phosim/
    #         notebooks/analysis_scripts/results_gaia/
    #         gMagGt11_w_2020_15/input


    # ensure that input/raw and input/rerun are empty 
    argString = '-rf  '+ os.path.join(outputDir,'input', 'rerun')
    print('Removing entire /input/rerun/* (-rf)')
    runProgram("rm", argstring=argString)

    argString = '-rf  '+ os.path.join(outputDir,'input', 'raw')
    print('Removing entire /input/raw/* (-rf)')
    runProgram("rm", argstring=argString)
    
    # remove files that are remade - this 
    # prevents problems with the registry ... 
    argString = os.path.join(outputDir, 'input')
    for remove in ['/isr*','/registry*' ,'/_mappe*']:
        print('Removing following files from input%s'%remove)
        runProgram("rm", argstring=argString+remove)
   
    # ensure that input/raw and input/rerun are empty 
    # print('Deleting content of input/raw/ and input/rerun/')
    # _eraseFolderContent(os.path.join(outputDir, 'input','raw'))
    # _eraseFolderContent(os.path.join(outputDir, 'input','rerun'))


    if defocalImg is True:
        print('We will make new defocal images in this run ')
        intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
        extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
        pertPath = os.path.join(outputDir, 'iter0','pert')
        if os.path.exists(intraPath):
            _eraseFolderContent(intraPath)
        if os.path.exists(extraPath):
            _eraseFolderContent(extraPath)
        if os.path.exists(pertPath):
            _eraseFolderContent(pertPath)




if not opd and not flats and not copy:
    print('Since copy=opd=flats=False, \
        we expect to see opd and calibs already there ... ') 

 # only do all that if not trying to just rerun the WFS ... 
# if not justWfs :  
#     # re-use the flats and calibs,
#     # since these have nothing to do with the actual stars ...
#     # -  copy from the sep_10  ...

#     # note : this is copying the results of previously ingested 
#     # calibration products. So to test a different version of the stack, 
#     # need to do the ingest ... 
#     if not opd and not flats and copy : 
#         print('Copying content of %s to re-use the flats and OPD files...'%copyDir)
        
#         # first copy the old results...
#         argString = '-a '+copyDir+'/. '+outputDir+'/'
#         runProgram("cp", argstring=argString)

#         # ensure that input/raw and input/rerun are empty 
#         print('Deleting content of input/raw/ and input/rerun/')
#         _eraseFolderContent(os.path.join(outputDir, 'input','raw'))
#         _eraseFolderContent(os.path.join(outputDir, 'input','rerun'))

#         # remove files that are remade
#         argString = os.path.join(outputDir, 'input')
#         runProgram("rm", argstring=argString+'/isr*')
#         runProgram("rm", argstring=argString+'/registry*')
#         runProgram("rm", argstring=argString+'/_mappe*')

#     # copy only opd and defocal images
#     if not opd and not defocalImg and copy : 
#         print('Copying content of %s to re-use the OPD images/files and defocal images ...'%copyDir)

#         # copy iter0/*
#         argString = '-a ' + copyDir+'/iter0/ ' + outputDir+'/iter0/'
#         iterDir = os.path.join(outputDir,'iter0')
#         if (not os.path.exists(iterDir)):
#             os.makedirs(iterDir)
#         runProgram("cp", argstring=argString)

#         # copy opdFiles 
#         argString = copyDir + '/opd* ' + outputDir+'/'
#         runProgram("cp", argstring=argString)



# Clobber
#if opd is True:
#    print('We will make new OPD files in this run')
#    _eraseFolderContent(outputDir)




# else:
#     if flats is True:
#         print('We will make new flats in this run')
#         _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
#         _eraseFolderContent(os.path.join(outputDir, 'input'))     
#     if defocalImg is True:
#         print('We will make new defocal images in this run ')
#         intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
#         extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
#         if os.path.exists(intraPath):
#             _eraseFolderContent(intraPath)
#         if os.path.exists(extraPath):
#             _eraseFolderContent(extraPath)


# read the star catalog from file ... 
# it conforms to the format expected by PhoSim 
# mv /data/epyc/users/suberlak/starCatGAIA.txt analysis_scripts/results_gaia
if selectSensors is 'comcam':
    skyFile = 'starCatGAIA_%s_%s.txt'%((field,catalogType))
    skyFilePath = os.path.join(topDir,skyFile)#starCatGAIA_gt11.txt')

# if selectSensors is 'wfs':
#     skyFilePath  = os.path.join(topDir, 'starCatGAIA_gt11_wfs.txt')
 
print('Using %s sensors and %s sky catalog'%(selectSensors,skyFilePath))
       
# set the opd.cmd and star.cmd files ...
opdCmd  = 'opdQuickBackground.cmd'
comcamCmd = 'starQuickBackground.cmd'

print('For PhoSim using /policy/cmdFile/%s and %s'%(opdCmd,comcamCmd))

# set the raInDeg,  decInDeg : 

# the center of field coords were first defined as 
# Galactic: 
from astropy.table import Table
from astropy.coordinates import SkyCoord
gt = Table(data=[['high','med','low','Baade'],
                            [0,0,0,1.02],
                           [85,40,10,-3.92 ]], 
                      names=['name', 'l_deg','b_deg'])

gaia_coords = SkyCoord(l=gt['l_deg'],b=gt['b_deg'], 
                       frame='galactic', unit='deg')
# convert them to equatorial 
gt['ra_deg']= gaia_coords.icrs.ra.deg
gt['dec_deg'] = gaia_coords.icrs.dec.deg

raInDeg = gt['ra_deg'][gt['name'] == field][0]
decInDeg = gt['dec_deg'][gt['name'] == field][0]
print('For this field, the raInDeg=%.3f, decInDeg=%.3f'%(raInDeg,decInDeg))


# initialize the baseComcamLoop.py Class 
print('Starting baseComcamLoop : ')
ccLoop = comcamLoop() 
ccLoop.main(phosimDir, numPro, iterNum, outputDir, testLabel, 
            isEimg=False,  genOpd=opd, genDefocalImg=defocalImg, genFocalImg=focalImg,
            genFlats=flats, useMinDofIdx=False,
            inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
            doDeblending=doDeblending, camDimOffset=camDimOffset, 
            postageImg=postageImg, opdCmdSettingsFile=opdCmd,
            comcamCmdSettingsFile=comcamCmd, selectSensors=selectSensors,
            splitWfsByMag =splitWfsByMag,raInDeg=raInDeg,decInDeg=decInDeg
            )
print('Done running ccLoop for GAIA \n\n')

# move the screenlog generated  by screen -LS  if it exists ... 
# it should be wherever the screen to run this code 
# got made ... 
screenlog_path = '/astro/store/epyc/users/suberlak/Commissioning/aos/ts_phosim/notebooks/analysis_scripts'
screenlog_default = os.path.join(screenlog_path,'screenlog.0')

if os.path.exists(screenlog_default):
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d_%H:%M:%S")
    screenlog = 'screenlog_%s_%s_%s.txt'%(topDir,expDir,date)
    argString = screenlog_default + ' '+outputDir+'/'+screenlog
    runProgram("mv ", argstring=argString)
    print('Screenlog saved as %s'%screenlog)
