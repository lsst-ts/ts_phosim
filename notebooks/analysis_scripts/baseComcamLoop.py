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

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt
from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, \
                                   getConfigDir
from lsst.ts.phosim.PlotUtil import plotFwhmOfIters


class baseComcamLoop():

    def main(self, phosimDir, numPro, iterNum, baseOutputDir,
            testName, isEimg=False, genOpd=True, genDefocalImg=True, 
            genFocalImg =True, genFlats=True,
            surveyFilter=None, starMag=15,
            useMinDofIdx=False, inputSkyFilePath="", m1m3ForceError=0.05,
            doDeblending=False, camDimOffset = None, postageImg=False,
            opdCmdSettingsFile='opdDefault.cmd',
            comcamCmdSettingsFile='starDefault.cmd', 
            instSettingFileName='starSingleExp.inst',
            selectSensors = 'comcam',
            splitWfsByMag=False, deblendDonutAlgo='convolveTemplate',
            centroidTemplateType='model', deblendTemplateType='isolatedDonutFromImage',
            raInDeg=None,decInDeg=None, rotAngInDeg=None):
        '''
        Code to run the full AOS loop on ComCam (or other sensors)

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
        comcamCmdSettingsFile: str, name of .cmd setting file for PhoSim when
            simulating the comcam images, should be located in
            /ts_phosim/policy/cmdFile/
        instSettingFileName : str, name of .inst setting file for PhoSim when 
            simulating images, should be located in /ts_phosim/policy/instFile/
        selectSensors: str, 'comcam'  for R22, or 'wfs' for corner wavefront sensors,
            a setting to pass explicitly to PhoSim  , also passed to _prepareOfcCalc,
            _prepareWepCalc
        splitWfsByMag: bool, whether to calculate the wfs for subsets of stars
            based on magnitude ranges, or not

        Parameters changing    ts_wep/policy/default.yaml  :
        --------------------------------------------------
        doDeblending : bool,  True by default
        camDimOffset : -150 , offset that is used to ignore stars that are that
            close to the CCD edge
        deblendDonutAlgo : str, a deblending algorithm to use if doDeblending=True,
            currently 'adapt' (old, pre-2020) or 'convolveTemplate' (new, May2020)
        templateType: str, which type of template to use with new centroid algorithms,
            'model', or 'phosim'

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

        useMinDofIdx: bool, True by default - whether to only use 10 hexapod
            positions and first 3 bending modes of M1M3 and M2


        '''
        
        # get the list of sensors  - by default it's comCam...
        sensorNameList = self._getComCamSensorNameList()

        # ... but it may be the wavefront sensing corner sensors ...
        if selectSensors is 'wfs':
            sensorNameList = self._getWfsSensorNameList()

        # Prepare the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            print('Making the calibration products ')
            # by default only make calibs for comcam

            if selectSensors is 'comcam':
                fakeFlatDir = self._makeCalibs(baseOutputDir, sensorNameList)

            if selectSensors is 'wfs':
                fakeFlatDir = self._makeCalibsWfs(baseOutputDir)

        # Make the ISR directory
        isrDirName = "input"
        isrDir = os.path.join(baseOutputDir, isrDirName)
        self._makeDir(isrDir)

        # Make the postage Image directory if needed
        if postageImg :
            postageImgDir  = os.path.join(baseOutputDir,'postage')
            self._makeDir(postageImgDir)
        else:
            postageImgDir = None

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
        print('Using the following settings for the telescope:')
        print('boresight (ra,dec) = %.3f,%.3f [deg]'%(raInDeg,decInDeg))
        print('rotation angle = %.3f [deg] '%rotAngInDeg)
        
        # Prepare the components
        print('Preparing the PhoSim component ')
        phosimCmpt = self._preparePhosimCmpt(phosimDir, filterType, raInDeg, decInDeg,
                                        rotAngInDeg, numPro, isEimg,
                                        m1m3ForceError)
        print('Preparing the wepCalc component ')
        wepCalc = self._prepareWepCalc(isrDir, filterType, raInDeg, decInDeg,
                                rotAngInDeg, isEimg, doDeblending, camDimOffset,
                                selectSensors,deblendDonutAlgo,centroidTemplateType,
                                deblendTemplateType)

        tele = phosimCmpt.getTele()
        defocalDisInMm = tele.getDefocalDistInMm()
        wepCalc.setDefocalDisInMm(defocalDisInMm)

        print('Preparing the ofcCalc component ')
        ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg,selectSensors)

        # Ingest the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            print('Ingesting calibration products')
            wepCalc.ingestCalibs(fakeFlatDir)

        # Only use 10 hexapod positions and first 3 bending modes of M1M3 and M2
        if (useMinDofIdx):
            self._useMinDofIdx(ofcCalc)

        # Set the telescope state to be the same as the OFC
        print('Setting the telescope state to be the same as OFC')
        state0 = ofcCalc.getStateAggregated()
        phosimCmpt.setDofInUm(state0)


        # decide which args should be added to PhoSim
        # they are prepended  - we add only working dir at this point
        # just prepend the working directory by default
        argPrepend = '-w ' + baseOutputDir+ ' '

        print('PhoSim added argPrepend is %s'%argPrepend)


        # Do the iteration
        obsId = 9006000
        opdZkFileName = str("opd.zer" + '.' + testName)
        wfsZkFileName = str("wfs.zer" + '.' + testName)
        opdPssnFileName = "PSSN.txt"
        outputDirName = "pert"
        outputImgDirName = "img"
        iterDefaultDirName = "iter"
        dofInUmFileName = "dofPertInNextIter.mat"
        skyInfoFileName = "skyComCamInfo.txt"
        for iterCount in range(iterNum):
            print('Starting iteration %d of %d'%(iterCount+1,iterNum))
            # Set the observation Id
            phosimCmpt.setSurveyParam(obsId=obsId)

            # The iteration directory
            iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

            # Set the output directory :   iter0/pert
            outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
            phosimCmpt.setOutputDir(outputDir)
            print('PhoSim outputDir is %s'%outputDir)

            # Set the output image directory:    iter0/img/
            outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                        outputImgDirName)
            phosimCmpt.setOutputImgDir(outputImgDir)
            print('PhoSim outputImgDir is %s'%outputImgDir)


            # Generate the OPD image
            if genOpd is True:
                t1 = datetime.datetime.now()

                if selectSensors is 'comcam':
                    argString = phosimCmpt.getComCamOpdArgsAndFilesForPhoSim(
                         cmdSettingFileName=opdCmdSettingsFile)

                # for LsstFamCam   == science sensors
                elif selectSensors is 'wfs':
                    argString = phosimCmpt.getLsstCamOpdArgsAndFilesForPhoSim(
                        cmdSettingFileName=opdCmdSettingsFile)

                argString = argPrepend + argString
                #argString = '-w $AOCLCOUTPUTPATH ' + argString
                print('Generating OPD with Phosim, argString is \n')
                print(argString)
                phosimCmpt.runPhoSim(argString)
                
                t2 = datetime.datetime.now()
                _print_duration(t2-t1)
            # Analyze the OPD data
            # this step creates iter0/img/PSSN.txt,
            # as well as opd.zer.xxx file
            # that describe the OPD
            print('Analyzing the OPD data ')

            if selectSensors is 'comcam':
                phosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
                                            pssnFileName=opdPssnFileName)
            elif selectSensors is 'wfs':
                phosimCmpt.analyzeLsstCamOpdData(zkFileName=opdZkFileName,
                                                pssnFileName=opdPssnFileName)

            # Get the PSSN from file
            pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
            print("   Calculated PSSN is %s." % pssn)

            # Get the GQ effective FWHM from file
            gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
            print("   GQ effective FWHM is %.4f." % gqEffFwhm)

            # Set the FWHM data
            listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
                                            opdPssnFileName, sensorNameList)
            ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)

            # Prepare the faked sky
            if (inputSkyFilePath == ""):
                # According to the OPD field positions
                metr = phosimCmpt.getOpdMetr()
                skySim = self._prepareSkySim(metr, starMag)
                print("Use the default OPD field positions to be star positions.")
                print("The star magnitude is chosen to be %.2f." % starMag)
            else:
                skySim = self._prepareSkySimBySkyFile(inputSkyFilePath)

            # Output the sky information 
            # xxx perhaps need to change that 
            # especially if using a star catalog ?? 
            skySim, wepCalc = self._outputSkyInfo(outputDir, skyInfoFileName,
                skySim, wepCalc)

            # Assign the entra- and intra-focal observation Id
            focalObsId = obsId
            extraObsId = obsId + 1
            intraObsId = obsId + 2
            

            #########################################
            # DEFOCAL IMAGES : GENERATE AND COLLECT
            ########################################


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
                if selectSensors is 'comcam':
                    rafts = ['22']
                    chips  = ['00','01','02',
                              '10','11','12',
                              '20','21','22']
                    sensors = ''
                    for r in rafts:
                        for c in chips:
                            s = "R%s_S%s|"%(r,c)
                            sensors += s
                    sensors = ' "%s" '%sensors

                if selectSensors is not None:
                    argPrepend +=  '-s ' + sensors+ ' '

                print('PhoSim added argPrepend is %s'%argPrepend)


                simSeed = 1000
                argStringList = phosimCmpt.getComCamStarArgsAndFilesForPhoSim(
                  extraObsId, intraObsId, skySim, simSeed=simSeed,
                  cmdSettingFileName=comcamCmdSettingsFile,
                  instSettingFileName=instSettingFileName)

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

                # Repackage the images based on the image type
                if (isEimg):
                    phosimCmpt.repackageComCamEimgFromPhoSim()
                else:
                    phosimCmpt.repackageComCamAmpImgFromPhoSim()


            
            # Collect the defocal images
            intraRawExpData = RawExpData()

            # it is iter0/img/intra/
            intraRawExpDir = os.path.join(outputImgDir,
                                        phosimCmpt.getIntraFocalDirName())
            intraRawExpData.append(intraObsId, 0, intraRawExpDir)

            extraRawExpData = RawExpData()
            # it is   iter0/img/extra/
            extraRawExpDir = os.path.join(outputImgDir,
                                        phosimCmpt.getExtraFocalDirName())
            extraRawExpData.append(extraObsId, 0, extraRawExpDir)


            #########################################
            # IN-FOCUS IMAGES : GENERATE AND COLLECT
            ########################################

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
                if selectSensors is 'comcam':
                    rafts = ['22']
                    chips  = ['00','01','02',
                              '10','11','12',
                              '20','21','22']
                    sensors = ''
                    for r in rafts:
                        for c in chips:
                            s = "R%s_S%s|"%(r,c)
                            sensors += s
                    sensors = ' "%s" '%sensors

                if selectSensors is not None:
                    argPrepend +=  '-s ' + sensors+ ' '

                print('\nPhoSim added argPrepend is %s'%argPrepend)


                simSeed = 1000
                argString = phosimCmpt.getComCamStarFocalPlaneArgsAndFilesForPhoSim(
                  obsId, skySim, simSeed=simSeed,
                  cmdSettingFileName=comcamCmdSettingsFile,
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
                else: 
                    phosimCmpt.repackageComCamAmpFocalImgFromPhoSim()

           
            # Collect the in-focus images
            focalRawExpData = RawExpData()

            # it is iter0/img/focal/
            focalRawExpDir = os.path.join(outputImgDir,
                                        phosimCmpt.getFocalDirName())
            focalRawExpData.append(focalObsId, 0, focalRawExpDir)


            ################################
            # CLEAR REGISTRY BEFORE INGEST 
            ################################

            # before ingesting ANY  images by WEP,  make sure that the previously ingested
            # ones are erased, especially in WFS-only mode !
            ingestedDir = os.path.join(isrDir, 'raw')
            if os.path.exists(ingestedDir):
                print('Removing the previously ingested raw images directory  %s \
                    before re-ingesting the images from iter0/img/...'%ingestedDir)
                argString = '-rf %s/'%ingestedDir
                runProgram("rm", argstring=argString)

            # also erase previously existing registry since this would mess the 
            # ingest process
            registryFile= os.path.join(isrDir,'registry.sqlite3')
            if os.path.exists(registryFile):
                print('Removing image registry file  %s '%registryFile)
                runProgram("rm", argstring=registryFile)


            ##################################
            # IN-FOCUS IMAGES : INGEST AND ISR   
            #################################

            # do the ingest and ISR on in-focus images : this is using 
            # just the beginning of     
            # wepCalc.calculateWavefrontErrors
         
            # When evaluating the eimage, the calibration products are not needed.
            # Therefore, need to make sure the camera mapper file exists.
            wepCalc._genCamMapperIfNeed()

            t1 = datetime.datetime.now()
            # Ingest the exposure data 
            print('Ingesting the in-focus images ')
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

            # Branch#1 : if we calculate wavefront errors 
            # for stars in magnitude bins,
            # first ingest and  do ISR on all images,
            # and then perform WFS calculation for each 
            # subset of target stars 
            if splitWfsByMag :
                print('Running WFS ingest and ISR once in split-stars-by-mag mode ')
                # an option to calculate WFS only
                # for stars of certain magnitude range  by
                # feeding the mag limits explicitly ...
                # Calculate the wavefront error and DOF

                ##############################
                #####   BEGIN PART 1    ######
                # calculateWavefrontErrors() in ts/wep/ctrlIntf/WEPCalculation.py , 
                # i.e. wepCalc.calculateWavefrontErrors()

                # first, ingest files just once ...
                rawExpData = intraRawExpData

                # When evaluating the eimage, the calibration products are not needed.
                # Therefore, need to make sure the camera mapper file exists.
                wepCalc._genCamMapperIfNeed()

                # Ingest the exposure data and do the ISR
                wepCalc._ingestImg(rawExpData)
                wepCalc._ingestImg(extraRawExpData)

                # Only the amplifier image needs to do the ISR
                imgType = wepCalc._getImageType()
                if (imgType == ImageType.Amp):
                    wepCalc._doIsr(isrConfigfileName="isr_config.py")

                # Set the butler inputs path to get the images
                butlerRootPath = wepCalc._getButlerRootPath()
                wepCalc.wepCntlr.setPostIsrCcdInputs(butlerRootPath)

                #####   END PART 1    ######
                ############################

                for lowMagnitude in [11,12,13,14,15]:
                    highMagnitude = lowMagnitude+1
                    print('    Calculating wavefront errors for stars between ')
                    print('    %d and %d magnitude'%(lowMagnitude,highMagnitude))

                    #############################
                    #####   BEGIN PART 2   ######
                    # wepCalc.calculateWavefrontErrors()
                    # Get the target stars map neighboring stars
                    neighborStarMap = wepCalc._getTargetStar(lowMagnitude=lowMagnitude,
                                                          highMagnitude=highMagnitude)

                    # Calculate the wavefront error
                    intraObsIdList = rawExpData.getVisit()
                    intraObsId = intraObsIdList[0]
                    if (extraRawExpData is None):
                        obsIdList = [intraObsId]
                    else:
                        extraObsIdList = extraRawExpData.getVisit()
                        extraObsId = extraObsIdList[0]
                        obsIdList = [intraObsId, extraObsId]

                    donutMap = wepCalc._calcWfErr(neighborStarMap, obsIdList,postageImg,postageImgDir)

                    listOfWfErr = wepCalc._populateListOfSensorWavefrontData(donutMap)

                    #####   END PART 2   ######
                    ###########################

                    ofcCalc.calculateCorrections(listOfWfErr)

                    zkFilenameAppend  = str(lowMagnitude)+'-'+str(highMagnitude)

                    # Record the wfs error with the same order as OPD for the comparison
                    phosimCmpt.reorderAndSaveWfErrFile(listOfWfErr, sensorNameList,
                                                   zkFileName=wfsZkFileName+zkFilenameAppend)



            # Branch #2 : if we want to calculate wavefront errors 
            # for all stars in the image, then use the 
            # built-in  wepCalc code, which performs 
            # ingest and  do ISR on all images, and 
            # then performs WFS calculation using all 
            # stars that fit selection criteria
            # 
            else: 
                print('Calculating the wavefront error ')
                  # Calculate the wavefront error and DOF
                t1 =datetime.datetime.now()
                _print_duration(t2-t1)

                listOfWfErr = wepCalc.calculateWavefrontErrors(
                    intraRawExpData, extraRawExpData=extraRawExpData,
                    postageImg=postageImg, postageImgDir = postageImgDir)
                ofcCalc.calculateCorrections(listOfWfErr)
                t2 =datetime.datetime.now()
                _print_duration(t2-t1)

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

        #sensorNameList = ["R22_S00", "R22_S01", "R22_S02", "R22_S10", "R22_S11",
        #                "R22_S12", "R22_S20", "R22_S21", "R22_S22"]

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

    def _getWfsSensorNameList(self):
        sensorNameList = ["R00_S22","R04_S20","R44_S00","R40_S02"]

        return sensorNameList

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
                        centroidTemplateType, deblendTemplateType):

        if (selectSensors is None) or (selectSensors is 'comcam'): # by default
            wepCalc = WEPCalculationFactory.getCalculator(CamType.ComCam, isrDirPath)
        elif selectSensors is 'wfs': # use LsstCam
            wepCalc = WEPCalculationFactory.getCalculator(CamType.LsstCam, isrDirPath)

        wepCalc.setFilter(filterType)
        wepCalc.setBoresight(raInDeg, decInDeg)
        wepCalc.setRotAng(rotAngInDeg)

        # call settingFile just once 
        settingFile = wepCalc.getSettingFile()
        
        # do all updates in order of appearance in policy/default.yaml
        if (isEimg):
            settingFile.updateSetting("imageType", "eimage")

        if bscDbType is not None:
            settingFile.updateSetting("bscDbType", bscDbType)

        if (doDeblending):
            settingFile.updateSetting("doDeblending", "True")
            settingFile.updateSetting("deblendDonutAlgo",deblendDonutAlgo)
            settingFile.updateSetting("centroidTemplateType", centroidTemplateType)
            settingFile.updateSetting("deblendTemplateType", deblendTemplateType)

        if camDimOffset  is not None :
            settingFile.updateSetting("camDimOffset", camDimOffset)
            
        # print info in order of appearance in policy/default.yaml
        print('Using following settings in ts_wep/policy/default.yaml:')
        print("imageType: %s"%settingFile.getSetting("imageType"))
        print("bscDbType: %s"%settingFile.getSetting("bscDbType"))
        print('camDimOffset: %s'% settingFile.getSetting("camDimOffset"))
        print("doDeblending:  %s"%settingFile.getSetting("doDeblending"))
        print("deblendDonutAlgo: %s"%settingFile.getSetting("deblendDonutAlgo"))
        print("centroidTemplateType: %s"%settingFile.getSetting("centroidTemplateType"))
        print("deblendTemplateType: %s"%settingFile.getSetting("deblendTemplateType"))

        return wepCalc


    def _prepareOfcCalc(self, filterType, rotAngInDeg, selectSensors):

        if (selectSensors is None) or (selectSensors is 'comcam'): # by default
            ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
        elif selectSensors is 'wfs': # use LsstCam
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


if __name__ == "__main__":

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop simulation (default is amp files).")
    parser.add_argument("--numOfProc", type=int, default=1,
                        help="number of processor to run PhoSim (default: 1)")
    parser.add_argument("--iterNum", type=int, default=5,
                        help="number of closed-loop iteration (default: 5)")
    parser.add_argument("--output", type=str, default="",
                        help="output directory")
    parser.add_argument("--testOutput", type=str, default="", help="test output directory")
    parser.add_argument("--testLabel", type=str, default="1",
                        help="filename identifier for test files")
    parser.add_argument('--eimage', default=False, action='store_true',
                        help='Use the eimage files')
    parser.add_argument('--minDof', default=False, action='store_true',
                        help='Use 10 hexapod positions and first 3 bending modes of M1M3 and M2')
    parser.add_argument("--skyFile", type=str, default="",
                        help="Star Id, ra, dec, and magnitude")
    parser.add_argument("--m1m3FErr", type=float, default=0.05,
                        help="Ratio of M1M3 actuator force error between 0 and 1 (default: 0.05)")
    parser.add_argument('--clobber', default=False, action='store_true',
                        help='Delete existing output directory')
    args = parser.parse_args()

    # Run the simulation
    phosimDir = getPhoSimPath()

    if (args.output == ""):
        outputDir = getAoclcOutputPath()
    else:
        outputDir = args.output

    os.makedirs(outputDir, exist_ok=True)
    if (args.clobber == True):
        _eraseFolderContent(outputDir)

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir

    ccLoop = baseComcamLoop()
    ccLoop.main(phosimDir, args.numOfProc, args.iterNum, outputDir, args.testLabel,
                isEimg=args.eimage, useMinDofIdx=args.minDof,
                inputSkyFilePath=args.skyFile, m1m3ForceError=args.m1m3FErr)
