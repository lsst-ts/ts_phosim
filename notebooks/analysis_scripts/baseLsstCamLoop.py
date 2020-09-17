#!/usr/bin/env python

import os
import argparse
import numpy as np
import shutil
import datetime

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import FilterType, CamType, runProgram
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
from lsst.utils import getPackageDir

class baseLsstCamLoop():

    def main(self, phosimDir, numPro, iterNum, baseOutputDir, 
            testName, isEimg=False, genOpd=True, genDefocalImg=True, genFlats=True,
            surveyFilter=None, 
            useMinDofIdx=False, inputSkyFilePath="", m1m3ForceError=0.05,
            doDeblending=False, camDimOffset = None, postageImg=False,
            opdCmdSettingsFile='opdDefault.cmd',
            starCmdSettingsFile='starDefault.cmd', 
            instSettingFile='starSingleExp.inst',
            selectSensors = 'wfs', 
            phosimRepackagerKeepOriginal=False,deblendDonutAlgo='convolveTemplate',
            centroidTemplateType='model', deblendTemplateType='isolatedDonutFromImage',
            bscDbType = 'file', db_filename  = 'bsc3.db3',
            raInDeg=None,decInDeg=None, rotAngInDeg=None):

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
        if db_filename  in os.listdir(bscDataDir):
            os.remove(os.path.join(bscDataDir,db_filename))
            print('\nRemoved old %s file'%db_filename)


        # get the list of sensors -  here the wavefront sensing corner sensors ... 
        if selectSensors is 'wfs':
            sensorNameList = self._getWfsSensorNameList()

        # Prepare the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            print('\nMaking calibration products ... ')
            #  make calibs for wfs 
            if selectSensors is 'wfs':
                fakeFlatDir = self._makeCalibsWfs(baseOutputDir)

        # Make the ISR directory
        isrDirName = "input"
        isrDir = os.path.join(baseOutputDir, isrDirName)
        if genFlats is True:
            self._makeDir(isrDir)

        # Make the postage Image directory if needed
        if postageImg :
            postageImgDir  = os.path.join(baseOutputDir,'postage')
            self._makeDir(postageImgDir)
        else:
            postageImgDir = None

        # Survey parameters - nothing comCam specific ...
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
        print('\nPreparing PhosimCmpt...')
        phosimCmpt = self._preparePhosimCmpt(phosimDir, filterType, raInDeg, decInDeg,
                                        rotAngInDeg, numPro, isEimg,
                                        m1m3ForceError)
        # here selectSensors =  'wfs'
        print('\nPreparing WepCalc...')
        wepCalc = self._prepareWepCalc(isrDir, filterType, raInDeg, decInDeg,
                                rotAngInDeg, isEimg, doDeblending, camDimOffset,
                                selectSensors,deblendDonutAlgo,centroidTemplateType,
                                deblendTemplateType)

        #tele = phosimCmpt.getTele()

        # NOTE : this step:
        
        # defocalDisInMm = tele.getDefocalDistInMm()
        # wepCalc.setDefocalDisInMm(defocalDisInMm)
        
        # is only possible for ComCam, where WepCalc 
        # is a nested instance of WEPCalculationOfComCam(WEPCalculationOfPiston),
        # WepCalculationOfPiston adds that method,
        # and WEPCalculationOfPiston(WEPCalculation)

        # WEPCalculationOfLsstCam(WEPCalculation),
        # so there is no such method 


        print('\nPreparing OfcCalc')
        ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg,selectSensors)

        # Ingest the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            print('\nIngesting calibration products ... ')
            wepCalc.ingestCalibs(fakeFlatDir)

        # Only use 10 hexapod positions and first 3 bending modes of M1M3 and M2
        if (useMinDofIdx):
            self._useMinDofIdx(ofcCalc)

        # Set the telescope state to be the same as the OFC
        state0 = ofcCalc.getStateAggregated()
        phosimCmpt.setDofInUm(state0)


        # Do the iteration
        obsId = 9006000
        opdZkFileName = str("opd.zer" + '.' + testName)
        wfsZkFileName = str("wfs.zer" + '.' + testName)
        opdPssnFileName = "PSSN.txt"
        outputDirName = "pert"
        outputImgDirName = "img"
        iterDefaultDirName = "iter"
        dofInUmFileName = "dofPertInNextIter.mat"
        skyInfoFileName = "skyInfo.txt"
        for iterCount in range(iterNum):

            # Set the observation Id
            phosimCmpt.setSurveyParam(obsId=obsId)

            # The iteration directory
            iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

            # Set the output directory
            outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
            phosimCmpt.setOutputDir(outputDir)
            print('\nPhoSim outputDir is %s'%outputDir)
            # Set the output image directory
            outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                        outputImgDirName)
            phosimCmpt.setOutputImgDir(outputImgDir)
            print('PhoSim outputImgDir is %s'%outputImgDir)

            # decide which args should be prepended to PhoSim 
            # just prepend the working directory by default 
            argPrepend = '-w ' + baseOutputDir+ ' ' 
            print('PhoSim added argPrepend is %s'%argPrepend)

            # Generate the OPD image
            # this makes iter0/img/opd_9006000_*.gz files,
            # each contains 255x255 array characterizing OPD 
            # in various locations (eg. for ComCam it was 
            # 1 per CCD, so 9 files,   for LSST Full Array Mode Cam - 
            # LsstFamCam, it is 31 locations in a ring, so 31 files )
            if genOpd is True:
                t1 = datetime.datetime.now()
                if selectSensors is 'wfs':
                    argString = phosimCmpt.getLsstFamCamOpdArgsAndFilesForPhoSim(
                            cmdSettingFileName=opdCmdSettingsFile)
                argString = argPrepend + argString

                print('\nGenerating OPD with Phosim, argString is \n')
                print(argString)
                phosimCmpt.runPhoSim(argString)
                t2 = datetime.datetime.now()
                _print_duration(t2-t1)

            # Analyze the OPD data

            # --> this takes the iter0/img/opd_**.gz files 

            # --> this makes the iter0/img/PSSN.txt  file that characterizes
            #     the OPD in point source sensitivity  - it contains 
            #     PSSN (one number per OPD file - for ComCam, it's 1 PSSN
            #     per CCD) , and Gaussian Quadrature  - derived FWHM (1 per 
            #     OPD file). 

            # --> this makes the opd.zer.xxx that characterizes OPD in 
            #     Zernikes  z4 to z22,  one row per input OPD file - for ComCam 
            #     that's 9 rows (one per CCD) ,  but for LsstFamCam - 31 rows 


            # if selectSensors is 'comcam':
            #     phosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
            #                                 pssnFileName=opdPssnFileName)
            if selectSensors is 'wfs':
                phosimCmpt.analyzeLsstFamCamOpdData(zkFileName=opdZkFileName,
                                            pssnFileName=opdPssnFileName)

            # Get the PSSN from file
            pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
            print("\nCalculated PSSN is %s." % pssn)

            # # Get the GQ effective FWHM from file
            gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
            print("GQ effective FWHM is %.4f." % gqEffFwhm)

            # Set the FWHM data
            listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
                                           opdPssnFileName, sensorNameList)
            ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)

            # Use the input sky catalog .... 
            skySim = self._prepareSkySimBySkyFile(inputSkyFilePath)

            # Output the sky information
            # this ensures that the star catalog
            # is always preserved in the directory with results 
            skySim, wepCalc = self._outputSkyInfo(outputDir, skyInfoFileName, 
                skySim, wepCalc)

            # Assign the entra- and intra-focal observation Id
            #extraObsId = obsId + 1
            intraObsId = obsId + 2

            # Generate the defocal images
            simSeed = 1000

            #if selectSensors is 'wfs':

            # There is actually nothing (that I can see) that is comcam - specific
            # here :  getComCamStarArgsAndFilesForPhoSim()
            # it gets self.tele.getDefocalDistInMm()  

            # We prepare OFC to be for LsstCam , 
            # and we ensure that phosim  has setDofInUm() as ofcCalc
            # thus I think that tele.getDefocalDistInMm() should have the 
            # WFS corner sensors values ... 
            
            print('\nUsing %s and %s for cmd and inst PhoSim files '%(starCmdSettingsFile,
                    instSettingFile))
            argString = phosimCmpt.getLsstCamStarArgsAndFilesForPhosim(
                 intraObsId=intraObsId, skySim=skySim, simSeed=simSeed,
                 cmdSettingFileName=starCmdSettingsFile,
                 instSettingFileName=instSettingFile)

            # Note, at this point PhosimCmpt.py  uses TeleFacade.py  
            # self.tele.writeStarInstFile() , which writes the 
            # camera configuration parameter 
            # camconfig  using the default 
            # self.sensorOn = {"sciSensorOn": True,
            #                  "wfSensorOn": True,
            #                  "guidSensorOn": False} 
            # i.e.  camconfig is 3, 
            # unless self.tele.sensorOn['sciSensorOn'] = False 
            # in which case camconfig is 2 ... 


            if genDefocalImg is True:

                # for running PhoSim on the defocal image on the selected WFS sensors
                # we prepend the sensors explicitly ... 

                # start again - prepend the working directory first 
                argPrepend = '-w ' + baseOutputDir+ ' ' 

                if selectSensors is 'wfs':
                    # note : need to select each half-chip ... 
                    sensors = "R00_S22_C0|R00_S22_C1|R04_S20_C0|R04_S20_C1|R44_S00_C0|R44_S00_C1|R40_S02_C0|R40_S02_C1"

                if selectSensors is not None: 
                    sensors = ' "%s" '%sensors # needed to pass the argument in comment signs 
                    argPrepend +=  '-s ' + sensors+ ' '

                print('\nFor generating defocal images with PhoSim, the argPrepend is ')
                print(argPrepend)

                argString = argPrepend + argString
                t1 = datetime.datetime.now()    
                print('Start time', t1.strftime("%Y-%m-%d_%H:%M:%S"))
                print('\n\nGenerating defocal images with Phosim, argString is ')
                print(argString)
                phosimCmpt.runPhoSim(argString)
                
                t2 = datetime.datetime.now()
                print('End time', t2.strftime("%Y-%m-%d_%H:%M:%S"))
                _print_duration(t2-t1)
                # Repackage the images based on the image type
                # Again, I don't see here anything comCam - specific...
                # the repackaging calls 
                # /epyc/projects/lsst_comm/phosim_utils/python/lsst/phosim/utils/phosim_repackager.py
            
                # just intra - focal : wrote a new function ...
                phosimCmpt.repackageLsstCamAmpImgFromPhosim(keepOriginal=phosimRepackagerKeepOriginal,
                                                            verbose=True) 

            # Collect the defocal images : only intra-focal...
          
            intraRawExpData = RawExpData() # instatiate the container class 
            intraRawExpDir = os.path.join(outputImgDir,
                                        phosimCmpt.getIntraFocalDirName())
            intraRawExpData.append(intraObsId, 0, intraRawExpDir)


            # before ingesting images by WEP,  make sure that the previously ingested 
            # ones are erased, especially in WFS-only mode !
            ingestedDir = os.path.join(isrDir, 'raw')
            if os.path.exists(ingestedDir):
                print('Removing the previously ingested raw images directory  %s \
                    before re-ingesting the images from iter0/img/...'%ingestedDir)
                argString = '-rf %s/'%ingestedDir
                runProgram("rm", argstring=argString)
            
            # also erase previously existing registry since this would mess the ingest process
            registryFile= os.path.join(isrDir,'registry.sqlite3')
            if os.path.exists(registryFile):
                print('Removing image registry file  %s '%registryFile)
                runProgram("rm", argstring=registryFile)


            # Calculate the wavefront error and DOF
            if selectSensors is 'wfs' : 
                sensorNameToIdFileName='sensorNameToIdWfs.yaml'
            else:
                sensorNameToIdFileName='sensorNameToId.yaml'
            print('Using sensor to ID translation from %s'%sensorNameToIdFileName)
            
            listOfWfErr = wepCalc.calculateWavefrontErrors(
                intraRawExpData, postageImg=postageImg, postageImgDir = postageImgDir,
                sensorNameToIdFileName=sensorNameToIdFileName)
            ofcCalc.calculateCorrections(listOfWfErr)

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

    def _getWfsSensorNameList(self):
        sensorNameList = ["R00_S22","R04_S20","R44_S00","R40_S02"]

        return sensorNameList

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

        if selectSensors is None: # by default
            ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
        elif selectSensors is 'wfs': # use LsstCam 
            ofcCalc = OFCCalculationFactory.getCalculator(InstName.LSST)

        ofcCalc.setFilter(filterType)
        ofcCalc.setRotAng(rotAngInDeg)
        ofcCalc.setGainByPSSN()

        return ofcCalc


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