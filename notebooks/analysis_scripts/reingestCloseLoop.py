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


class baseReingest():

    def main(self, baseOutputDir, genFlats=True,  isEimg=False, ingestDefocal=True, ingestFocal=True,
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

        # Prepare the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            print('Making the calibration products ')
            # by default only make calibs for comcam
            fakeFlatDir = self._makeCalibs(baseOutputDir, sensorNameList)


        # Make the ISR directory ( it was cleaned before )
        isrDirName = "input"
        isrDir = os.path.join(baseOutputDir, isrDirName)
        self._makeDir(isrDir)

     
        # Survey parameters, taken from ts_phosim/policy/surveySettings.yaml
        surveySettingFilePath = os.path.join(getConfigDir(),
                                            "surveySettings.yaml")
        surveySettings = ParamReader(filePath=surveySettingFilePath)
  
        filterType = FilterType.fromString(
                surveySettings.getSetting("filterType"))
       
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
        # print('Preparing the PhoSim component ')
        # phosimCmpt = self._preparePhosimCmpt(phosimDir, filterType, raInDeg, decInDeg,
        #                                 rotAngInDeg, numPro, isEimg,
        #                                 m1m3ForceError)
        print('Preparing the wepCalc component ')
        wepCalc = self._prepareWepCalc(isrDir, filterType, raInDeg, decInDeg,
                                rotAngInDeg, isEimg)

        # tele = phosimCmpt.getTele()
        # defocalDisInMm = tele.getDefocalDistInMm()
        # wepCalc.setDefocalDisInMm(defocalDisInMm)

   
        # Ingest the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            print('Ingesting calibration products')
            wepCalc.ingestCalibs(fakeFlatDir)

     


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

        if ingestDefocal : 
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
            wepCalc._ingestImg(intraRawExpData)
            wepCalc._ingestImg(extraRawExpData)

            # Only the amplifier image needs to do the ISR
            imgType = wepCalc._getImageType()
            if (imgType == ImageType.Amp):
                print('Performing ISR on defocal images ')
                wepCalc._doIsr(isrConfigfileName="isr_config.py")



        #########################################
        # IN-FOCUS IMAGES :  COLLECT, INGEST AND ISR  
        ########################################

        if ingestFocal : 
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
            wepCalc._ingestImg(focalRawExpData)
        

        # Only the amplifier image needs to do the ISR
        # but we're only doing amplifier images for 
        # in-focus images ... 
        if isEimg:
            print("No need to do the ISR on in-focus e-images ")
            pass 
        else: 
            print('Performing the ISR on in-focus amp images ')
            wepCalc._doIsr(isrConfigfileName="isr_config.py")

    
          



               

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
                        isEimg):

        wepCalc = WEPCalculationFactory.getCalculator(CamType.ComCam, isrDirPath)

        wepCalc.setFilter(filterType)
        wepCalc.setBoresight(raInDeg, decInDeg)
        wepCalc.setRotAng(rotAngInDeg)

        if (isEimg):
            settingFile = wepCalc.getSettingFile()
            settingFile.updateSetting("imageType", "eimage")

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

