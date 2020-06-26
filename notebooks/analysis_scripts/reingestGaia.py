import os
import datetime
import numpy as np
from reingestCloseLoop import baseReingest  as baseReingest
from reingestCloseLoop import _eraseFolderContent
from createPhosimCatalogNew import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir
from lsst.ts.wep.Utility import runProgram
from astropy.table import Table
from astropy.coordinates import SkyCoord



# origin directory  
field = 'high'
inputDir = 'results_gaia/dr2_%s_gt11'%field

# a directory in which we store the rerun
outputDir = 'results_gaia/dr2_%s_gt11_reingest'%field


if not os.path.exists(outputDir):
    print('Copying entire content of %s ...'%inputDir)
    print('This includes calibs, OPD, defocal images, ingest, etc.')

    #first copy the old results...
    argString = '-a '+inputDir+'/. '+outputDir+'/'
    runProgram("cp", argstring=argString)


    # ensure that input/ is empty
    argString = '-rf  '+ os.path.join(outputDir,'input')
    print('Removing entire /input/* (-rf)')
    runProgram("rm", argstring=argString)


print('\nStarting reingest  ')

print('The outputDir is %s'%outputDir)


# set the raInDeg,  decInDeg : 

# the center of field coords were first defined as 
# Galactic: 

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


# initialize the reingestCloseLoop.py Class 
print('Starting reingestCloseLoop : ')
reingest = baseReingest() 
reingest.main(baseOutputDir = outputDir, genFlats=True,  isEimg=False, 
                  ingestDefocal=True, ingestFocal=True,
                  raInDeg=raInDeg,decInDeg=decInDeg
                  )
print('Done running reingest for %s \n\n'%outputDir)

# move the screenlog generated  by screen -LS  if it exists ... 
screenlog_path = '/astro/store/epyc/users/suberlak/Commissioning/aos/ts_phosim/notebooks/analysis_scripts'
screenlog_default = os.path.join(screenlog_path,'screenlog.0')

if os.path.exists(screenlog_default):
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d_%H:%M:%S")
    screenlog = 'screenlog_%s.txt'%(date)
    argString = screenlog_default + ' '+outputDir+'/'+screenlog
    runProgram("mv ", argstring=argString)
    print('Screenlog saved as %s'%screenlog)
