#!/usr/bin/env python

import os
import argparse
import numpy as np
import shutil
import datetime

from lsst.ts.wep.Utility import FilterType, CamType, runProgram
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory
from lsst.ts.wep.ctrlIntf.RawExpData import RawExpData


class baseReingest():

    def main(self, baseOutputDir, genFlats=True, ingestRawDefocal=True, isrDefocal = True,
             ingestRawFocal=True, isrFocal = True, isrConfigfileName = "isr_config.py",
             rerunName = 'run1'):
        '''
        Code to redo the ingestion process for the AOS loop products. This assumes that eg. 
        there was an AOS loop run with one version of the stack (eg. 2020_15),
        and we want to remake the calibs, ingest them to make butler calibRegistry  and input/flats,
        ingest the raw images to input/raw,  making butler registry , 
        and perform the ISR, making input/rerun/run1 ...
        
        Another scenario is if we first run the AOS loop with certain ISR settings, 
        and want to rerun the ISR with different settings. With the update of 
        makeGainImages.py, which makes fake_flats, it is recommended 
        to remake the calibs, ingest them, and then re-do the ISR. 


        Parameters:
        ------------
     
        baseOutputDir: str, the base output directory path for all output, eg.
            'results_gaia/gMagGt11_w_2020_15_test'
    

        genFlats: boolean,  True/False  - whether to generate with PhoSim the
            calibration files, in /fake_flats/

        ingestRawDefocal : boolean,  True/False  - whether to ingest the raw defocal images
        isrDefocal :  boolean,  True/False  - whether to redo the ISR on raw defocal images 
      
        ingestRawFocal : boolean,  True/False  - whether to ingest the raw in-focus images
        isrFocal :  boolean,  True/False  - whether to redo the ISR on raw in-focus images 

        isrConfigfileName : str, 'isr_config.py' by default : filename where the ISR 
            settings should be saved 

        rerunName : str, 'run1' is the default, this kwarg is passed to wep_calc._doIsr(), and 
            supersedes the settingFile.getSetting("rerunName") - if we want to keep the 
            original postISR image products,  use eg. 'run2' 
    
        '''
        
        # get the list of sensors  - by default it's comCam...
        sensorNameList = self._getComCamSensorNameList()

        # Prepare the calibration products (only for the amplifier images)
        if genFlats :
            print('Making the calibration products ')
            # by default only make calibs for comcam
            t1 = datetime.datetime.now()
            fakeFlatDir = self._makeCalibs(baseOutputDir, sensorNameList)
            t2 = datetime.datetime.now() 
            _print_duration(t2-t1)

        # Make the ISR directory ( if needed )
        isrDirName = "input"
        isrDir = os.path.join(baseOutputDir, isrDirName)
        self._makeDir(isrDir)

        # initialize wep calc      
        print('Preparing the wepCalc component ')
        wepCalc = WEPCalculationFactory.getCalculator(CamType.ComCam, isrDir)
   
        # Ingest the calibration products (only for the amplifier images)
        if genFlats  :
            print('Ingesting calibration products')
            t1 = datetime.datetime.now() 
            wepCalc.ingestCalibs(fakeFlatDir)
            t2 = datetime.datetime.now() 
            _print_duration(t2-t1)

        ####################
        #     ITERATION 
        ####################

        obsId = 9006000

        outputImgDirName = "img"
        outputDirName = "pert"
        iterDefaultDirName = "iter"
        iterCount = 0 
            

        # The iteration directory
        iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

        # Set the output image directory:    iter0/img/
        outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                    outputImgDirName)



        # Assign the entra- and intra-focal observation Id
        focalObsId = obsId
        extraObsId = obsId + 1
        intraObsId = obsId + 2
        

        #########################################
        # DEFOCAL IMAGES: COLLECT, INGEST AND ISR 
        ########################################

        if ingestRawDefocal : 
            print('Collecting defocal images')

            # Collect the defocal images
            intraRawExpData = RawExpData()

            # it is iter0/img/intra/
            intraRawExpDir = os.path.join(outputImgDir,'intra')
            intraRawExpData.append(intraObsId, 0, intraRawExpDir)

            extraRawExpData = RawExpData()
            # it is   iter0/img/extra/
            extraRawExpDir = os.path.join(outputImgDir,'extra')
            extraRawExpData.append(extraObsId, 0, extraRawExpDir)

            # When evaluating the eimage, the calibration products are not needed.
            # Therefore, need to make sure the camera mapper file exists.
            wepCalc._genCamMapperIfNeed()

            print('Ingesting defocal images')
            # Ingest the exposure data and do the ISR
            t1 = datetime.datetime.now() 
            wepCalc._ingestImg(intraRawExpData)
            wepCalc._ingestImg(extraRawExpData)
            t2= datetime.datetime.now() 
            _print_duration(t2-t1)

        if isrDefocal:
            # Only the amplifier image needs to do the ISR
            #imgType = wepCalc._getImageType()
            #if (imgType == ImageType.Amp):
            print('Performing ISR on defocal amplifier images ')
            t1 = datetime.datetime.now()
            wepCalc._doIsr(isrConfigfileName=isrConfigfileName,
                           rerunName = rerunName)
            t2 = datetime.datetime.now()
            _print_duration(t2-t1)



        #########################################
        # IN-FOCUS IMAGES :  COLLECT, INGEST AND ISR  
        ########################################

        if ingestRawFocal : 
            # Collect the in-focus images
            focalRawExpData = RawExpData()

            # it is iter0/img/focal/
            focalRawExpDir = os.path.join(outputImgDir, 'focal')
            focalRawExpData.append(focalObsId, 0, focalRawExpDir)

            # When evaluating the eimage, the calibration products are not needed.
            # Therefore, need to make sure the camera mapper file exists.
            wepCalc._genCamMapperIfNeed()

            # Ingest the exposure data 
            print('Ingesting the in-focus images ')
            t1 = datetime.datetime.now() 
            wepCalc._ingestImg(focalRawExpData)
            t2 = datetime.datetime.now() 
            _print_duration(t2-t1)
        
        if isrFocal:
            # Only the amplifier image needs to do the ISR
            # but we're only doing amplifier images for 
            # in-focus images ... 
            # if isEimg:
            #     print("No need to do the ISR on in-focus e-images ")
            #     pass 
            # else: 
            print('Performing the ISR on in-focus amp images ')
            t1 = datetime.datetime.now() 
            wepCalc._doIsr(isrConfigfileName=isrConfigfileName,
                           rerunName = rerunName)
            t2 = datetime.datetime.now() 
            _print_duration(t2-t1)


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

