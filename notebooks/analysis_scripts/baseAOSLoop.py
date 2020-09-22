#!/usr/bin/env python

import os
import argparse
import numpy as np
import shutil
import datetime

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import FilterType, CamType, runProgram, ImageType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory
from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData
from lsst.ts.wep.ctrlIntf.MapSensorNameAndId import MapSensorNameAndId

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, \
								   getConfigDir
from lsst.ts.phosim.PlotUtil import plotFwhmOfIters
from lsst.utils import getPackageDir


class baseAOSLoop():

	def main(self, phosimDir, numPro, iterNum, baseOutputDir,
			testName, isEimg=False, genOpd=True, genDefocalImg=True, 
			genFocalImg =True, genFlats=True,
			surveyFilter=None, starMag=15,
			useMinDofIdx=False, inputSkyFilePath="", m1m3ForceError=0.05,
			doDeblending=False, camDimOffset = None, postageImg=False,
			opdCmdSettingsFile='opdDefault.cmd',
			starCmdSettingsFile='starDefault.cmd', 
			instSettingFileName='starSingleExp.inst',
			selectSensors = 'comcam',
			deblendDonutAlgo='convolveTemplate',
			centroidTemplateType='model', deblendTemplateType='isolatedDonutFromImage',
			bscDbType = 'file', dbFileName  = 'bsc2.db3', expWcs = False,
			raInDeg=None,decInDeg=None, rotAngInDeg=None):
		'''
		Code to run the full AOS loop on ComCam / LsstCam / LsstFamCam

		Parameters:
		------------
		phosimDir : str, a directory with phosim.py,  returned by
			lsst.ts.phosim.Utility.getPhoSimPath . For UW, it is
			'/epyc/projects/lsst_comm/phosim_syseng4/'
		numPro : int , 10 by default - number of processors used for
			parallel calculation by PhoSim
			NB. - since each raytrace is itself massively parallel, using eg. 10
			cores may use 35 in effect (%use often shows eg. 350% per CPU...)

		baseOutputDir: str, the base output directory path for all output, eg.
			'results_gaia/gMagGt11_w_2020_15_test'
		testName: str , a label for the test appended to wfs.zer , opd.zer  files,
			 eg 'gaia' yields  `wfs.zer.gaia',  'opd.zer.gaia'....

		genOpd: boolean,  True/False  - whether to generate the Optical
			Path Difference files, in  /iter0/img/opd_9006000_*.fits.gz
		genDefocalImg: boolean,  True/False  - whether to generate with PhoSim
			the defocal images, in /iter0/img/extra/  and  iter0/img/intra/
		genFocalImg: boolean, True/False  - whether to generate with PhoSim
			the in-focus images, in /iter0/img/focal/   
		genFlats: boolean,  True/False  - whether to generate with PhoSim the
			calibration files, in /fake_flats/

		inputSkyFilePath: str, path to the input star catalog (with
			ID | RA | DEC |  MAG  ),  eg.  '/results_gaia/starCatalog.txt'
		postageImg: bool, True by default  - whether to save postage stamp
			images of stars during the ts_wep  calcuation of wavefront error.
			They are saved in eg.  /results_gaia/gMagGt11/postage/

		opdCmdSettingsFile: str, name of .cmd setting file for PhoSim when
			simulating the OPD images, should be located in
			/ts_phosim/policy/cmdFile/
		starCmdSettingsFile: str, name of .cmd setting file for PhoSim when
			simulating the comcam images, should be located in
			/ts_phosim/policy/cmdFile/
		instSettingFileName : str, name of .inst setting file for PhoSim when 
			simulating images, should be located in /ts_phosim/policy/instFile/
		selectSensors: str, 'comcam'  for R22, or 'wfs' for corner wavefront sensors,
			a setting to pass explicitly to PhoSim  , also passed to _prepareOfcCalc,
			_prepareWepCalc

		Parameters changing    ts_wep/policy/default.yaml  :
		--------------------------------------------------
		doDeblending : bool,  True by default
		camDimOffset : -150 , offset that is used to ignore stars that are that
			close to the CCD edge
		deblendDonutAlgo : str, a deblending algorithm to use if doDeblending=True,
			currently 'adapt' (old, pre-2020) or 'convolveTemplate' (new, May2020)
		templateType: str, which type of template to use with new centroid algorithms,
			'model', or 'phosim'
		expWcs: str, False by default (using WcsSol), if True - using PhosimWcsSol
		dbFileName: str, bsc_db1.db3 by default , name of database file to use 
		bscDbType: str, 'file' by default ('refCat', 'image' are other options)

		Parameters changing ts_phosim/policy/surveySettings.yaml : 
		--------------------------------------------------------
		surveyFilter: str, by default: None, which means that we read the value from
			surveySettings.yaml for "filterType : ref", which means reference filter,
			and in this context its LSST g-filter 
		raInDeg: float, by default: None, which means that 0.0 is read from setting file 
		decInDeg: float, by default: None, which means that 0.0 is read from setting file
		rotAngInDeg: float, by default: None, which means that 0.0 is read from setting file


		Parameters not often changed (legacy):
		--------------------------------------
		starMag: int, a magnitude of a test star, if there isn't a catalog of
			sources provided
		iterNum: int, number of iterations - 1 by default ... I've never changed that
		m1m3ForceError: int, 0.05 by default (why?)
		isEimg: bool, False by default - whether to make an electronic or amplifier
				image.

		useMinDofIdx: bool, False by default - whether to only use 10 hexapod
			positions and first 3 bending modes of M1M3 and M2


		'''
		
		# NB : the bscDbType has to be 
		# updated before wep_calc gets initialized.
		# Once wep_calc is called , there is no way 
		# to change bscDbType (any .updateSetting("bscDbType") call 
		# will not take desired effect)
		path_to_ts_wep = getPackageDir("ts_wep")
		setting_filename = 'default.yaml'
		path_to_setting_file = os.path.join(path_to_ts_wep, 'policy',setting_filename)
		settingFile = ParamReader(filePath=path_to_setting_file)
		bscDbTypeInFile = settingFile.getSetting("bscDbType")
		print('\n%s contains : '%path_to_setting_file, 
			  bscDbTypeInFile)
		if bscDbTypeInFile != bscDbType:
			# In the following we update the setting for bscDbType,
			# saving the change in the default.yaml file 
			settingFile.updateSetting("bscDbType", bscDbType)
			settingFile.saveSetting(filePath=path_to_setting_file)

			# check that the change indeed x`took place 
			settingFile = ParamReader(filePath=path_to_setting_file)
			print('After change: ', settingFile.getSetting("bscDbType"))
		
		bscDataDir = os.path.join(path_to_ts_wep, 'tests/testData')
		if dbFileName  in os.listdir(bscDataDir):
			os.remove(os.path.join(bscDataDir,dbFileName))
			print('\nRemoved old %s file'%dbFileName)


		# update expWcs setting ... 
		settingFile= ParamReader(filePath=path_to_setting_file)
		settingFile.updateSetting("expWcs", expWcs)
		settingFile.saveSetting(filePath=path_to_setting_file)
		settingFile = ParamReader(filePath=path_to_setting_file)
		print('After change: expWcs : ', settingFile.getSetting("expWcs"))


		# Prepare the calibration products (only for the amplifier images)
		if ((not isEimg) & (genFlats is True)):
			print('\nMaking the calibration products ')
			# by default only make calibs for comcam

			if selectSensors is 'comcam':
				sensorNameList = self._getComCamSensorNameList()
				fakeFlatDir = self._makeCalibs(baseOutputDir, sensorNameList)

			elif selectSensors is 'lsstfamcam':
				sensorNameList = self._getLsstFamCamSensorNameList()
				fakeFlatDir = self._makeCalibs(baseOutputDir, sensorNameList)

			elif selectSensors is 'lsstcam':
				fakeFlatDir = self._makeCalibsWfs(baseOutputDir)


		# Make the ISR directory
		isrDirName = "input"
		isrDir = os.path.join(baseOutputDir, isrDirName)
		self._makeDir(isrDir)

		
		# Survey parameters, taken from ts_phosim/policy/surveySettings.yaml
		surveySettingFilePath = os.path.join(getConfigDir(),
											"surveySettings.yaml")
		surveySettings = ParamReader(filePath=surveySettingFilePath)
		if surveyFilter is None:
			filterType = FilterType.fromString(
				surveySettings.getSetting("filterType"))
		else:
			filterType = FilterType.fromString(surveyFilter)
		if raInDeg is None:
			raInDeg = surveySettings.getSetting("raInDeg")
		if decInDeg is None:
			decInDeg = surveySettings.getSetting("decInDeg")
		if rotAngInDeg is None :
			rotAngInDeg = surveySettings.getSetting("rotAngInDeg")
		print('\nUsing the following settings for the telescope:')
		print('boresight (ra,dec) = %.3f,%.3f [deg]'%(raInDeg,decInDeg))
		print('rotation angle = %.3f [deg] '%rotAngInDeg)
		
		# Prepare the components
		print('\nPreparing the PhoSim component ')
		phosimCmpt = self._preparePhosimCmpt(phosimDir, filterType, raInDeg, decInDeg,
										rotAngInDeg, numPro, isEimg,
										m1m3ForceError)

		print('\nPreparing the wepCalc component ')
		wepCalc = self._prepareWepCalc(isrDir, filterType, raInDeg, decInDeg,
								rotAngInDeg, isEimg, doDeblending, camDimOffset,
								selectSensors,deblendDonutAlgo,centroidTemplateType,
								deblendTemplateType, dbFileName)
		
		# this step is not possible for the WFS sensors 
		# since there the defocal images are achieved by the inherent property
		# of 1500 micron z-axis offset  for the corner sensors, as recorder in 
		# focalplanelayout.txt file in phosim_syseng4/data/lsst/ and ts_wep/policy/ 

		if (selectSensor == 'comcam' ) or (selectSensor == 'lsstfamcam' ):
			tele = phosimCmpt.getTele()
			defocalDisInMm = tele.getDefocalDistInMm()
			wepCalc.setDefocalDisInMm(defocalDisInMm)

		print('\nPreparing the ofcCalc component ')
		ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg,selectSensors)

		# Ingest the calibration products (only for the amplifier images)
		if ((not isEimg) & (genFlats is True)):
			print('\nIngesting calibration products')
			wepCalc.ingestCalibs(fakeFlatDir)

		# Only use 10 hexapod positions and first 3 bending modes of M1M3 and M2
		if (useMinDofIdx):
			self._useMinDofIdx(ofcCalc)

		# Set the telescope state to be the same as the OFC
		print('\nSetting the telescope state to be the same as OFC')
		state0 = ofcCalc.getStateAggregated()
		phosimCmpt.setDofInUm(state0)


		# Do the iteration
		obsId = 9006000
		opdZkFileName = str("opd.zer" + '.' + testName)
		wfsZkFileName = str("wfs.zer" + '.' + testName)
		opdPssnFileName = "PSSN.txt"
		outputDirName = "pert"
		outputImgDirName = "img"
		outputPostageDirName = "postage"
		iterDefaultDirName = "iter"
		dofInUmFileName = "dofPertInNextIter.mat"
		skyInfoFileName = "skyInfo.txt"

		for iterCount in range(iterNum):
			print('\nStarting iteration %d of %d'%(iterCount+1,iterNum))
			# Set the observation Id
			phosimCmpt.setSurveyParam(obsId=obsId)

			# The iteration directory
			iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

			# Set the output directory :   iter0/pert
			outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
			phosimCmpt.setOutputDir(outputDir)
			print('\nPhoSim outputDir is %s'%outputDir)

			# Set the output image directory:    iter0/img/
			outputImgDir = os.path.join(baseOutputDir, iterDirName,
										outputImgDirName)
			phosimCmpt.setOutputImgDir(outputImgDir)
			print('\nPhoSim outputImgDir is %s'%outputImgDir)

			# Make the postage image directory if needed:
			# new for each iteration ... 
			if postageImg :
				outputPostageDir  = os.path.join(baseOutputDir,iterDirName,
										outputPostageDirName)
				self._makeDir(outputPostageDir)
			else:
				outputPostageDir = None


			# decide which args should be added to PhoSim
			# they are prepended  - we add only working dir at this point
			# just prepend the working directory by default
			argPrepend = '-w ' + baseOutputDir+ ' '

			print('\nPhoSim added argPrepend is %s'%argPrepend) 

			# Generate the OPD image
			if genOpd is True:
				t1 = datetime.datetime.now()
				print('Start time', t1.strftime("%Y-%m-%d_%H:%M:%S"))
				
				# generate OPD at 9 locations 
				if selectSensors == 'comcam':
					print('\nGenerating OPD at 9 field locations corresponding to\
						centers of ComCam sensors.')
					argString = phosimCmpt.getComCamOpdArgsAndFilesForPhoSim(
						 cmdSettingFileName=opdCmdSettingsFile)

				# generate OPD at 31 field locations 
				if (selectSensors == 'lsstfamcam' ) or (selectSensors == 'lsstcam'):
					print('\nGenerating OPD at 31 field locations spread across the\
						full focal plane.')
					argString = phosimCmpt.getLsstFamCamOpdArgsAndFilesForPhoSim(
						cmdSettingFileName=opdCmdSettingsFile)

				argString = argPrepend + argString
				print('Generating OPD with Phosim, argString is \n')
				print(argString)
				phosimCmpt.runPhoSim(argString)
				
				t2 = datetime.datetime.now()
				print('End time', t2.strftime("%Y-%m-%d_%H:%M:%S"))
				_print_duration(t2-t1)

			# Analyze the OPD data
			# this step creates iter0/img/PSSN.txt,
			# as well as opd.zer.xxx file
			# that describe the OPD

			# --> this takes the iter0/img/opd_**.gz files 

			# --> this makes the iter0/img/PSSN.txt  file that characterizes
			#     the OPD in point source sensitivity  - it contains 
			#     PSSN (one number per OPD file - for ComCam, it's 1 PSSN
			#     per CCD) , and Gaussian Quadrature  - derived FWHM (1 per 
			#     OPD file). 

			# --> this makes the opd.zer.xxx that characterizes OPD in 
			#     Zernikes  z4 to z22,  one row per input OPD file - for ComCam 
			#     that's 9 rows (one per CCD) ,  but for LsstFamCam - 31 rows 

			print('\nAnalyzing the OPD data ')

			if selectSensors == 'comcam':
				phosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
											pssnFileName=opdPssnFileName)

			elif (selectSensors == 'lsstcam') or (selectSensors == 'lsstfamcam'):
				phosimCmpt.analyzeLsstFamCamOpdData(zkFileName=opdZkFileName,
												pssnFileName=opdPssnFileName)

			# Get the PSSN from file
			pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
			print("   Calculated PSSN is %s." % pssn)

			# Get the GQ effective FWHM from file
			gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
			print("   GQ effective FWHM is %.4f." % gqEffFwhm)

			# Set the FWHM data providing the list of sensors on which OPD was 
			# evaluated 
			if selectSensors  == 'comcam':
				sensorNameList = self._getComCamSensorNameList()

			if (selectSensors == 'lsstcam' ) or  (selectSensors == 'lsstfamcam'):  
				# need to choose 31 of 189 sensors for which OPD was evaluated 
				# we use for that the sensorNameToFieldIdx.yaml
				# translation file from ts_ofc/policy/lsst/
				sensorNameList = self._getUniqueLsstFamCamOpdSensors()

			listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
											opdPssnFileName, sensorNameList)
			ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)

			# Prepare the faked sky
			if (inputSkyFilePath == ""):
				# According to the OPD field positions
				metr = phosimCmpt.getOpdMetr()
				skySim = self._prepareSkySim(metr, starMag)
				print("\nUse the default OPD field positions to be star positions.")
				print("The star magnitude is chosen to be %.2f." % starMag)
			else:
				skySim = self._prepareSkySimBySkyFile(inputSkyFilePath)

			# Output the sky information  - this is used if 
			# bscDbType :  file 
			skySim, wepCalc = self._outputSkyInfo(outputDir, skyInfoFileName,
				skySim, wepCalc)

			# Assign the entra- and intra-focal observation Id
			focalObsId = obsId
			extraObsId = obsId + 1
			intraObsId = obsId + 2
			

			#########################################
			# DEFOCAL IMAGES : GENERATE AND COLLECT
			########################################

			simSeed = 1000

			# Generate the defocal images
			if genDefocalImg is True:


				# Ensure that the folder is empty 
				# especially if copying some files 
				intraRawExpDir = os.path.join(outputImgDir,
										phosimCmpt.getIntraFocalDirName())
				if os.path.exists(intraRawExpDir): # iter0/img/extra/
					print('Before proceeding, cleaned up %s '%intraRawExpDir)
					_eraseFolderContent(intraRawExpDir)
				
				extraRawExpDir = os.path.join(outputImgDir,
										phosimCmpt.getExtraFocalDirName())
				if os.path.exists(extraRawExpDir): # iter0/img/extra/
					print('Before proceeding, cleaned up %s '%extraRawExpDir)
					_eraseFolderContent(extraRawExpDir)


				# just prepend the working directory by default
				argPrepend = '-w ' + baseOutputDir+ ' '

				# then prepend argument to run PhoSim only on R22
				# for defocal images
				if selectSensors == 'comcam':
					sensorNameList = self._getComCamSensorNameList()

				elif selectSensors == 'lsstcam':
					sensorNameList = self._getLsstCamSensorNameList()

				elif selectSensors == 'lsstfamcam':
					sensorNameList = self._getLsstFamCamSensorNameList()

				
				sensorNameString =  self._sensorNameListToString(sensorNameList)

				if selectSensors is not None:
					argPrepend +=  '-s  "%s"  '%sensorNameString

				print('\n PhoSim added argPrepend is %s'%argPrepend)


				# for comcam or lsstfamcam it's the same 
				if (selectSensors == 'comcam' ) or (selectSensors == 'lsstfamcam') :
					argStringList = phosimCmpt.getComCamStarArgsAndFilesForPhoSim(
							  extraObsId, intraObsId, skySim, simSeed=simSeed,
							  cmdSettingFileName=starCmdSettingsFile,
							  instSettingFileName=instSettingFileName)

				# corner sensor only create 1 "set" of images 
				elif selectSensors == 'lsstcam':
					argString = phosimCmpt.getLsstCamStarArgsAndFilesForPhosim(
							  intraObsId=intraObsId, skySim=skySim, simSeed=simSeed,
							  cmdSettingFileName=starCmdSettingsFile,
							  instSettingFileName=instSettingFileName)
					argStringList  = [argString]


				# iterate over extra and intra-focal images 
				for argString in argStringList:
					t1 = datetime.datetime.now()    
					print('Start time', t1.strftime("%Y-%m-%d_%H:%M:%S"))
					argString = argPrepend + argString
					print('Generating defocal images with Phosim\n')
					print(argString)
					phosimCmpt.runPhoSim(argString)

					t2 = datetime.datetime.now()
					print('End time', t2.strftime("%Y-%m-%d_%H:%M:%S"))
					_print_duration(t2-t1)
				
				if (selectSensors == 'comcam') or (selectSensors == 'lsstfamcam') :
					# Repackage the images based on the image type
					if (isEimg):
						phosimCmpt.repackageComCamEimgFromPhoSim()
					else:
						phosimCmpt.repackageComCamAmpImgFromPhoSim()
				

				elif selectSensors == 'lsstcam' :
					# just intra - focal : wrote a new function ...
					phosimCmpt.repackageLsstCamAmpImgFromPhosim(
						keepOriginal=phosimRepackagerKeepOriginal,
						verbose=True) 

			if (selectSensors == 'comcam') or (selectSensors == 'lsstfamcam') :
				# Collect the defocal images
				intraRawExpData = RawExpData()

				# it is iter0/img/intra/
				intraRawExpDir = os.path.join(outputImgDir, phosimCmpt.getIntraFocalDirName())
				intraRawExpData.append(intraObsId, 0, intraRawExpDir)

				extraRawExpData = RawExpData()
				# it is   iter0/img/extra/
				extraRawExpDir = os.path.join(outputImgDir, phosimCmpt.getExtraFocalDirName())
				extraRawExpData.append(extraObsId, 0, extraRawExpDir)
			 
			elif selectSensors == 'lsstcam' : # only intra-focal imgs to collect 
				intraRawExpData = RawExpData() 
				intraRawExpDir = os.path.join(outputImgDir, phosimCmpt.getIntraFocalDirName())
				intraRawExpData.append(intraObsId, 0, intraRawExpDir)
	 
			#########################################
			# IN-FOCUS IMAGES : GENERATE AND COLLECT
			########################################

			if (selectSensors == 'comcam') or (selectSensors == 'lsstfamcam') :
				if genFocalImg is True : 

					# Ensure that the folder is empty 
					# especially if copying some files 
					focalRawExpDir = os.path.join(outputImgDir,
											phosimCmpt.getFocalDirName())
					if os.path.exists(focalRawExpDir): # iter0/img/focal/
						print('Before proceeding, cleaned up %s '%focalRawExpDir)
						_eraseFolderContent(focalRawExpDir)


				   # just prepend the working directory by default
					argPrepend = '-w ' + baseOutputDir+ ' '


					# then prepend argument to run PhoSim only on R22
					# just in case some stars provided to PhoSim
					# had streaks or extended 

					if selectSensors == 'comcam':
						sensorNameList = self._getComCamSensorNameList()

					elif selectSensors == 'lsstfamcam':
						sensorNameList = self._getLsstFamCamSensorNameList()

					sensorNameString =  self._sensorNameListToString(sensorNameList)

					argPrepend +=  '-s  "%s"  '%sensorNameString

					print('\nPhoSim added argPrepend is %s'%argPrepend)

					argString = phosimCmpt.getComCamStarFocalPlaneArgsAndFilesForPhoSim(
					  obsId, skySim, simSeed=simSeed,
					  cmdSettingFileName=starCmdSettingsFile,
					  instSettingFileName=instSettingFileName)

					argString = argPrepend + argString
					print('\nGenerating focal plane images with Phosim\n')
					
					t1 = datetime.datetime.now()
					print('Start time', t1.strftime("%Y-%m-%d_%H:%M:%S"))
					print(argString)
					
					phosimCmpt.runPhoSim(argString)

					t2 = datetime.datetime.now()
					print('End time', t2.strftime("%Y-%m-%d_%H:%M:%S"))
					_print_duration(t2-t1)

					# Repackage the images : these are amp images 
					# so  I only make a  function for amp images 
					if (isEimg):
						print("Repackaging function for in-focus e-images \
							   doesn't exist yet")
						pass 
					else:  # its the same for lsstfamcam ... 
						phosimCmpt.repackageComCamAmpFocalImgFromPhoSim()

			   
					# Collect the in-focus images
					focalRawExpData = RawExpData()

					# it is iter0/img/focal/
					focalRawExpDir = os.path.join(outputImgDir,
												phosimCmpt.getFocalDirName())
					focalRawExpData.append(focalObsId, 0, focalRawExpDir)



			##################################
			# IN-FOCUS IMAGES : INGEST AND ISR   
			#################################
			if (selectSensors == 'comcam') or (selectSensors == 'lsstfamcam') :
				if genFocalImg is True :
					# do the ingest and ISR on in-focus images : this is using 
					# just the beginning of     
					# wepCalc.calculateWavefrontErrors
				 
					# When evaluating the eimage, the calibration products are not needed.
					# Therefore, need to make sure the camera mapper file exists.
					wepCalc._genCamMapperIfNeed()

					t1 = datetime.datetime.now()
					# Ingest the exposure data 
					print('\nIngesting the in-focus images ')
					wepCalc._ingestImg(focalRawExpData)
					t2 =datetime.datetime.now()
					_print_duration(t2-t1)

					# Only the amplifier image needs to do the ISR
					# but we're only doing amplifier images for 
					# in-focus images ... 
					if isEimg:
						print("No need to do the ISR on in-focus e-images ")
						pass 
					else: 
						print('Performing the ISR on in-focus amp images ')
					   wepCalc._doIsr(isrConfigfileName="isr_config.py")

		
			########################################
			# DEFOCAL IMAGES : INGEST, ISR, WEPCALC
			########################################

			print('\nCalculating the wavefront error ')
			t1 =datetime.datetime.now()

			listOfWfErr = wepCalc.calculateWavefrontErrors(
				intraRawExpData, extraRawExpData=extraRawExpData,
				postageImg=postageImg, postageImgDir = outputPostageDir)

			ofcCalc.calculateCorrections(listOfWfErr)
			t2 =datetime.datetime.now()
			_print_duration(t2-t1)


			if selectSensors == 'comcam':
				sensorNameList = self._getComCamSensorNameList()

			elif selectSensors  == 'lsstcam':
				sensorNameList  = ["R00_S22", "R04_S20","R40_S02","R44_S00"]

			elif selectSensors == 'lsstfamcam':
				sensorNameList == self._getLsstFamCamSensorNameList()


			# Record the wfs error with the same order as OPD for the comparison
			phosimCmpt.reorderAndSaveWfErrFile(listOfWfErr, sensorNameList,
										   zkFileName=wfsZkFileName)

		   
			# Set the new aggregated DOF to phosimCmpt
			dofInUm = ofcCalc.getStateAggregated()
			phosimCmpt.setDofInUm(dofInUm)

			# Save the DOF file
			phosimCmpt.saveDofInUmFileForNextIter(
				dofInUm, dofInUmFileName=dofInUmFileName)

			# Add the observation ID by 10 for the next iteration
			obsId += 10

		# Summarize the FWHM
		pssnFiles = [os.path.join(baseOutputDir, "%s%d" % (iterDefaultDirName, num),
					outputImgDirName, opdPssnFileName) for num in range(iterNum)]
		saveToFilePath = os.path.join(baseOutputDir, "fwhmIters.png")
		plotFwhmOfIters(pssnFiles, saveToFilePath=saveToFilePath)

	def _outputSkyInfo(self, outputDir, skyInfoFileName, skySim, wepCalc):

		outputSkyInfoFilePath = os.path.join(outputDir, skyInfoFileName)
		skySim.exportSkyToFile(outputSkyInfoFilePath)
		wepCalc.setSkyFile(outputSkyInfoFilePath)

		return skySim, wepCalc

	def _getComCamSensorNameList(self):

		chips  = ['00','01','02',
			  '10','11','12',
			  '20','21','22']
		rafts = ['22']
		sensors = []
		for r in rafts:
			for c in chips:
				s = "R%s_S%s"%(r,c)
				sensors.append(s)
		sensorNameList = sensors
		return sensorNameList
	
	def _getLsstFamCamSensorNameList(self):
		# I assume it includes the ComCam 
		chips  = ['00','01','02',
				  '10','11','12',
				  '20','21','22']
		rafts = ['14','24','34', 
			'03','13','23','33','43',
			'02','12','22','32','42',
			'01','11','21','31','41',
				 '10','20','30']
		sensors = []
		for r in rafts:
			for c in chips:
				s = "R%s_S%s"%(r,c)
				sensors.append(s)
		sensorNameList = sensors
		return sensorNameList

	def _getLsstCamSensorNameList(self):
		sensorNameList = ["R00_S22_C0","R00_S22_C1", 
						  "R04_S20_C0","R04_S20_C1",
						  "R40_S02_C0","R40_S02_C1",
						  "R44_S00_C0","R44_S00_C0"
						 ]

		return sensorNameList

	def _sensorNameListToString(self,sensorNameList):
		sensors  = ''
		for sensor in sensorNameList:
			sensors += "%s|"%sensor
		return sensors

	def _makeCalibs(self, outputDir, sensorNameList):

		fakeFlatDirName = "fake_flats"
		fakeFlatDir = os.path.join(outputDir, fakeFlatDirName)
		self._makeDir(fakeFlatDir)

		detector = " ".join(sensorNameList)
		self._genFakeFlat(fakeFlatDir, detector)

		return fakeFlatDir

	def _makeCalibsWfs(self,outputDir):
		fakeFlatDirName = "fake_flats"
		fakeFlatDir = os.path.join(outputDir, fakeFlatDirName)
		self._makeDir(fakeFlatDir)

		currWorkDir = os.getcwd()
		os.chdir(fakeFlatDir)
		command = "makeGainImages.py"
		argstring = "--just_wfs"
		runProgram(command,argstring=argstring)
		os.chdir(currWorkDir)

		return fakeFlatDir

	def _makeDir(self, directory):

		if (not os.path.exists(directory)):
			os.makedirs(directory)


	def _genFakeFlat(self, fakeFlatDir, detector):

		currWorkDir = os.getcwd()

		os.chdir(fakeFlatDir)
		self._makeFakeFlat(detector)
		os.chdir(currWorkDir)


	def _makeFakeFlat(self, detector):

		command = "makeGainImages.py"
		argstring = "--detector_list %s" % detector
		runProgram(command, argstring=argstring)


	def _preparePhosimCmpt(self, phosimDir, filterType, raInDeg, decInDeg, rotAngInDeg,
						numPro, isEimg, m1m3ForceError):

		# Set the Telescope facade class
		tele = TeleFacade()
		tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
		tele.setPhoSimDir(phosimDir)
		print('Using phosim.py located in %s'%phosimDir)

		# Prepare the phosim component
		phosimCmpt = PhosimCmpt(tele)

		# Set the telescope survey parameters
		boresight = (raInDeg, decInDeg)
		zAngleInDeg = 27.0912
		phosimCmpt.setSurveyParam(filterType=filterType, boresight=boresight,
								zAngleInDeg=zAngleInDeg, rotAngInDeg=rotAngInDeg)

		# Update the setting file if needed
		settingFile = phosimCmpt.getSettingFile()
		if (numPro > 1):
			settingFile.updateSetting("numPro", numPro)
		if isEimg:
			settingFile.updateSetting("e2ADC", 0)

		# Set the seed number for M1M3 surface
		seedNum = 6
		phosimCmpt.setSeedNum(seedNum)

		# Set the M1M3 force error
		phosimCmpt.setM1M3ForceError(m1m3ForceError)

		return phosimCmpt


	def _prepareWepCalc(self, isrDirPath, filterType, raInDeg, decInDeg, rotAngInDeg,
						isEimg,doDeblending, camDimOffset, selectSensors,deblendDonutAlgo,
						centroidTemplateType, deblendTemplateType, dbFileName):

		if (selectSensors is None) or (selectSensors == 'comcam'): # by default
			wepCalc = WEPCalculationFactory.getCalculator(CamType.ComCam, isrDirPath)

		elif selectSensors is 'lsstcam': # use LsstCam
			wepCalc = WEPCalculationFactory.getCalculator(CamType.LsstCam, isrDirPath)

		elif selectSensors is 'lsstfamcam':
			wepCalc = WEPCalculationFactory.getCalculator(CamType.LsstFamCam, isrDirPath)

		wepCalc.setFilter(filterType)
		wepCalc.setBoresight(raInDeg, decInDeg)
		wepCalc.setRotAng(rotAngInDeg)

		# call settingFile just once 
		settingFile = wepCalc.getSettingFile()
		
		# do all updates in order of appearance in policy/default.yaml
		if (isEimg):
			settingFile.updateSetting("imageType", "eimage")

		if (doDeblending):
			settingFile.updateSetting("doDeblending", "True")
			settingFile.updateSetting("deblendDonutAlgo",deblendDonutAlgo)
			settingFile.updateSetting("centroidTemplateType", centroidTemplateType)
			settingFile.updateSetting("deblendTemplateType", deblendTemplateType)
		dbRelativePath = 'tests/testData/%s'%dbFileName
		settingFile.updateSetting("defaultBscPath", dbRelativePath)

		if camDimOffset  is not None :
			settingFile.updateSetting("camDimOffset", camDimOffset)
			
		# print info in order of appearance in policy/default.yaml
		print('\nUsing following settings in ts_wep/policy/default.yaml:')
		print("imageType: %s"%settingFile.getSetting("imageType"))
		print("bscDbType: %s"%settingFile.getSetting("bscDbType"))
		print('camDimOffset: %s'% settingFile.getSetting("camDimOffset"))
		print("doDeblending:  %s"%settingFile.getSetting("doDeblending"))
		print("deblendDonutAlgo: %s"%settingFile.getSetting("deblendDonutAlgo"))
		print("centroidTemplateType: %s"%settingFile.getSetting("centroidTemplateType"))
		print("deblendTemplateType: %s"%settingFile.getSetting("deblendTemplateType"))
		print("defaultBscPath: %s"%settingFile.getSetting("defaultBscPath"))

		return wepCalc


	def _prepareOfcCalc(self, filterType, rotAngInDeg, selectSensors):

		if (selectSensors is None) or (selectSensors == 'comcam'): # by default
			ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)

		elif (selectSensors == 'lsstcam') or (selectSensors == 'lsstfamcam'):
			ofcCalc = OFCCalculationFactory.getCalculator(InstName.LSST)

		ofcCalc.setFilter(filterType)
		ofcCalc.setRotAng(rotAngInDeg)
		ofcCalc.setGainByPSSN()

		return ofcCalc


	def _prepareSkySim(self, opdMetr, starMag):

		skySim = SkySim()

		starId = 0
		raInDegList, declInDegList = opdMetr.getFieldXY()
		for raInDeg, declInDeg in zip(raInDegList, declInDegList):
			# It is noted that the field position might be < 0. But it is not the
			# same case for ra (0 <= ra <= 360).
			if (raInDeg < 0):
				raInDeg += 360.0
			skySim.addStarByRaDecInDeg(starId, raInDeg, declInDeg, starMag)
			starId += 1

		return skySim


	def _prepareSkySimBySkyFile(self, inputSkyFilePath):

		skySim = SkySim()

		absSkyFilePath = os.path.abspath(inputSkyFilePath)
		skySim.addStarByFile(absSkyFilePath)

		return skySim


	def _useMinDofIdx(self, ofcCalc):

		ztaac = ofcCalc.getZtaac()

		m1m3Bend = np.zeros(20, dtype=int)
		m1m3Bend[0: 3] = 1

		m2Bend = np.zeros(20, dtype=int)
		m2Bend[0: 3] = 1

		ztaac.setZkAndDofInGroups(m1m3Bend=m1m3Bend, m2Bend=m2Bend)

	def _eraseFolderContent(targetDir):

		for theFile in os.listdir(targetDir):
			filePath = os.path.join(targetDir, theFile)
			if os.path.isfile(filePath):
				os.unlink(filePath)
			elif os.path.isdir(filePath):
				shutil.rmtree(filePath)


	def _print_duration(delta):
		''' Convenience function to print execution time 
		between time1 and time2. 
		
		Parameters:
		----------
		delta : datetime.timedelta() object,
		result of eg. 
		time1 = datetime.datetime.now()
		time2 = datetime.datetime.now()
		delta = time2 - time1 
		
		Returns:
		--------
		None

		'''
		delta_sec = delta.total_seconds()
		delta_min = delta_sec / 60
		delta_hr = delta_min / 60
		print('    It took %.3f minutes, i.e. %.5f hours ' % (delta_min, delta_hr))

	def _getUniqueLsstFamCamOpdSensors(self):
		refSensorNameList = self._getLsstFamCamSensorNameList()
		sensorNameToIdFileName  = 'sensorNameToId.yaml'
		mapSensorNameAndId = MapSensorNameAndId(sensorNameToIdFileName)

		path_to_ts_ofc = getPackageDir("ts_ofc")
		mappingFilePath = os.path.join(path_to_ts_ofc , 'policy/lsst', 'sensorNameToFieldIdx.yaml')
		_mappingFile = ParamReader()
		_mappingFile.setFilePath(mappingFilePath)

		fieldIdx = []
		for sensor in refSensorNameList:
			field = _mappingFile.getSetting(sensor)
			fieldIdx.append(int(field))
			
		uniqueFieldIdx = np.unique(fieldIdx)
		uniqueFieldIdxLt31 = uniqueFieldIdx[uniqueFieldIdx<31]

		#This shows all sensors corresponding to each fieldIdx 
		oneSensorPerFieldIdx = []
		for field in uniqueFieldIdxLt31:
			#print(field, np.array(refSensorNameList)[fieldIdx == field])
			oneSensorPerFieldIdx.append(np.array(refSensorNameList)[fieldIdx == field][0])


if __name__ == "__main__":

	# Set the parser
	parser = argparse.ArgumentParser(
		description="Run AOS closed-loop simulation (default is amp files).")
	parser.add_argument("--numPro", type=int, default=1,
						help="number of processors to run PhoSim (default: 1)")
	parser.add_argument("--iterNum", type=int, default=5,
						help="number of closed-loop iteration (default: 5)")
	parser.add_argument("--outputDir", type=str, default="",
						help="output directory")
	parser.add_argument("--testLabel", type=str, default="1",
						help="filename identifier for test files")
	parser.add_argument('--isEimg', default=False, action='store_true',
						help='Use the eimage files')
	parser.add_argument('--genFocalImg', default=False, action='store_true',
						help='Generate in-focus images')
	parser.add_argument('--useMinDofIdx', default=False, action='store_true',
						help='Use 10 hexapod positions and first 3 bending modes of M1M3 and M2')
	parser.add_argument("--skyFile", type=str, default="",
						help="Star Id, ra, dec, and magnitude")
	parser.add_argument("--m1m3ForceError", type=float, default=0.05,
						help="Ratio of M1M3 actuator force error between 0 and 1 (default: 0.05)")
	parser.add_argument('--clobber', default=False, action='store_true',
						help='Delete existing output directory')
	args = parser.parse_args()

	# Run the simulation
	phosimDir = getPhoSimPath()

	os.makedirs(args.outputDir, exist_ok=True)
	if (args.clobber == True):
		_eraseFolderContent(args.outputDir)

	ccLoop = baseComcamLoop()
	ccLoop.main(phosimDir, args.numPro, args.iterNum, args.outputDir, args.testLabel,
				isEimg=args.isEimg, genFocalImg = args.genFocalImg, 
				useMinDofIdx=args.useMinDofIdx,
				inputSkyFilePath=args.skyFile, m1m3ForceError=args.m1m3ForceError)
	
