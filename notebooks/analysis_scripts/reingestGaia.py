import os
import datetime
import numpy as np
import argparse
from astropy.table import Table
from astropy.coordinates import SkyCoord
from reingestCloseLoop import baseReingest  as baseReingest
from lsst.ts.wep.Utility import runProgram


parser = argparse.ArgumentParser(
    description='Reingest results of running ComCam AOS closed loop for GAIA DR2 catalogs')
parser.add_argument("--field", type=str, default="",
                    help='name of GAIA field, eg. high, med, low, Baade. If other, leave empty. ')

parser.add_argument("--inputDir", type=str, default="",
                    help='input directory' )
parser.add_argument("--outputDir", type=str, default="",
                    help='output directory' )
parser.add_argument("--raInDeg", type=float, default=0, 
                    help='RA of the field center in degrees')
parser.add_argument("--decInDeg", type=float, default=0, 
                    help='DEC of the field center in degrees')

args = parser.parse_args()



# Try setting the origin and destination 
# directories from the field name 
field = args.field 
# if len(field)>0:
#     print('Based on provided field name %s, we are assuming the following:'%field)
#     # Set origin directory 
#     inputDir = 'results_gaia/dr2_%s_gt11'%field

#     # a directory in which we store the rerun
#     outputDir = 'results_gaia/dr2_%s_gt11_reingest'%field
# else:

if (len(args.inputDir)>0) and (len(args.outputDir)>0) : 
    print("\nUsing the following as the inputDir and outputDir")
    inputDir = 'results_gaia/%s'%args.inputDir
    outputDir = 'results_gaia/%s'%args.outputDir
    print(inputDir)
    print(outputDir)
     
else:
	print("Need to provide inputDir and outputDir names relative to results_gaia/")


if len(field)>0:
    print('Based on provided field name %s, we are assuming the following:'%field)
    
    # set the raInDeg,  decInDeg ,
    # if the field name is provided 
    # the center of field coords were first defined as 
    # Galactic: 
    print('\nCalculating raInDeg, decInDeg based on field name ')
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

elif abs(args.raInDeg - 0.0) > 1e-3  :
	    print("Reading the raInDeg, decInDeg for the field center from parser")
	    raInDeg = args.raInDeg
	    decInDeg =args.decInDeg 

else:

    raise ValueError('Need to provide a named field, or the input raInDeg and declInDeg ')

# Copy the results of the prior PhoSim run
if not os.path.exists(outputDir):
    print('\nCopying entire content of %s ...'%inputDir)
    print('This includes calibs, OPD, defocal images, \
        ingested images, postISR, etc.')

    #first copy the old results...
    argString = '-a '+inputDir+'/. '+outputDir+'/'
    runProgram("cp", argstring=argString)


    # ensure that input/ is empty
    argString = '-rf  '+ os.path.join(outputDir,'input')
    print('Removing entire /input/* (-rf)')
    runProgram("rm", argstring=argString)


print('\nStarting reingest  ')

print('The outputDir is %s'%outputDir)



# initialize the reingestCloseLoop.py Class 
print('Starting reingestCloseLoop : ')
reingest = baseReingest() 
reingest.main(baseOutputDir = outputDir, genFlats=True,  isEimg=False, 
                  ingestDefocal=True, ingestFocal=True,
                  raInDeg=raInDeg,decInDeg=decInDeg
                  )
print('Done running reingest for %s \n\n'%outputDir)

# # move the screenlog generated  by screen -LS  if it exists ... 
screenlog_path = '/astro/store/epyc/users/suberlak/Commissioning/aos/ts_phosim/notebooks/analysis_scripts'
screenlog_default = os.path.join(screenlog_path,'screenlog.0')

if os.path.exists(screenlog_default):
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d_%H:%M:%S")
    screenlog = 'screenlog_%s.txt'%(date)
    argString = screenlog_default + ' '+outputDir+'/'+screenlog
    runProgram("mv ", argstring=argString)
    print('Screenlog saved as %s'%screenlog)
