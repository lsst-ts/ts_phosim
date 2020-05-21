#!/usr/bin/env python

import os
import argparse
import numpy as np
import shutil

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


class AnalysisPhosimCmpt(PhosimCmpt):

    def _writeOpdZkFile(self, zkFileName, rotOpdInDeg):
        """Write the OPD in zk file.
        OPD: optical path difference.
        Parameters
        ----------
        zkFileName : str
            OPD in zk file name.
        rotOpdInDeg : float
            Rotate OPD in degree in the counter-clockwise direction.
        """

        testOutputDir = os.environ["closeLoopTestDir"]
        filePath = os.path.join(testOutputDir, zkFileName)
        opdData = self._mapOpdToZk(rotOpdInDeg)
        header = "The followings are OPD in rotation angle of %.2f degree in um from z4 to z22:" % (
            rotOpdInDeg)
        np.savetxt(filePath, opdData, header=header)

    def reorderAndSaveWfErrFile(self, listOfWfErr, refSensorNameList,
                                zkFileName="wfs.zer"):
        """Reorder the wavefront error in the wavefront error list according to
        the reference sensor name list and save to a file.
        The unexisted wavefront error will be a numpy zero array. The unit is
        um.
        Parameters
        ----------
        listOfWfErr : list [lsst.ts.wep.ctrlIntf.SensorWavefrontData]
            List of SensorWavefrontData object.
        refSensorNameList : list
            Reference sensor name list.
        zkFileName : str, optional
            Wavefront error file name. (the default is "wfs.zer".)
        """

        # Get the sensor name that in the wavefront error map
        wfErrMap = self._transListOfWfErrToMap(listOfWfErr)
        nameListInWfErrMap = list(wfErrMap.keys())

        # Reorder the wavefront error map based on the reference sensor name
        # list.
        reorderedWfErrMap = dict()
        for sensorName in refSensorNameList:
            if sensorName in nameListInWfErrMap:
                wfErr = wfErrMap[sensorName]
            else:
                numOfZk = self.getNumOfZk()
                wfErr = np.zeros(numOfZk)
            reorderedWfErrMap[sensorName] = wfErr

        # Save the file
        testOutputDir = os.environ["closeLoopTestDir"]
        filePath = os.path.join(testOutputDir, zkFileName)
        wfsData = self._getWfErrValuesAndStackToMatrix(reorderedWfErrMap)
        header = "The followings are ZK in um from z4 to z22:"
        np.savetxt(filePath, wfsData, header=header)

class baseWfsWep():

    def main(self, phosimDir, numPro, iterNum, baseOutputDir, 
            testName, isEimg=False, genOpd=True, genDefocalImg=True, genFlats=True,
            surveyFilter=None, starMag=15, 
            useMinDofIdx=False, inputSkyFilePath="", m1m3ForceError=0.05,
            doDeblending=False, camDimOffset = None, postageImg=False,
            opdCmdSettingsFile='opdDefault.cmd',
            comcamCmdSettingsFile='starDefault.cmd', selectSensors = 'wfs'):

        # get the list of sensors  - by default it's comCam...
        # sensorNameList = self._getComCamSensorNameList()

        # # ... but it may be the wavefront sensing corner sensors ... 
        if selectSensors is 'wfs':
            sensorNameList = self._getWfsSensorNameList()

        # Prepare the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            #  make calibs for wfs 
            # if selectSensors is 'comcam':
            #     fakeFlatDir = self._makeCalibs(baseOutputDir, sensorNameList)
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
        raInDeg = surveySettings.getSetting("raInDeg")
        decInDeg = surveySettings.getSetting("decInDeg")
        rotAngInDeg = surveySettings.getSetting("rotAngInDeg")

        # Prepare the components
        print('\nPreparing PhosimCmpt...')
        phosimCmpt = self._preparePhosimCmpt(phosimDir, filterType, raInDeg, decInDeg,
                                        rotAngInDeg, numPro, isEimg,
                                        m1m3ForceError)
        # here selectSensors =  'wfs'
        print('\nPreparing WepCalc...')
        wepCalc = self._prepareWepCalc(isrDir, filterType, raInDeg, decInDeg,
                                rotAngInDeg, isEimg, doDeblending, camDimOffset,
                                selectSensors)

        tele = phosimCmpt.getTele()

        # this step only possible for ComCam, where WepCalc 
        # is a nested instance of WEPCalculationOfComCam(WEPCalculationOfPiston),
        # WepCalculationOfPiston adds that method,
        # and WEPCalculationOfPiston(WEPCalculation)

        # WEPCalculationOfLsstCam(WEPCalculation),
        # so there is no such method 

        #defocalDisInMm = tele.getDefocalDistInMm()
        #wepCalc.setDefocalDisInMm(defocalDisInMm)
        print('\nPreparing OfcCalc')
        ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg,selectSensors)

        # Ingest the calibration products (only for the amplifier images)
        if ((not isEimg) & (genFlats is True)):
            wepCalc.ingestCalibs(fakeFlatDir)

        # Only use 10 hexapod positions and first 3 bending modes of M1M3 and M2
        if (useMinDofIdx):
            self._useMinDofIdx(ofcCalc)

        # Set the telescope state to be the same as the OFC
        state0 = ofcCalc.getStateAggregated()
        phosimCmpt.setDofInUm(state0)


        # decide which args should be added to PhoSim 
        # they are prepended 
        # this applies both to OPD and to star image 

        # just prepend the working directory by default 
        argPrepend = '-w ' + baseOutputDir+ ' ' 
       

        # then prepend argument to run PhoSim only on R22 
        # if selectSensors is 'comcam':  
        #     rafts = ['22']
        #     chips  = ['00','01','02', 
        #               '10','11','12',
        #               '20','21','22']
        #     sensors = ''
        #     for r in rafts:
        #         for c in chips:
        #             s = "R%s_S%s|"%(r,c) 
        #             sensors += s 

        if selectSensors is 'wfs':
            sensors = "R00_S22|R04_S20|R44_S00|R40_S02"
            #rafts = ['00','04', '40', '44']
            #chips  = ['00','01','02', 
            #          '10','11','12',
            #          '20','21','22']
            

        if selectSensors is not None: 
            sensors = ' "%s" '%sensors # needed to pass the argument in 
            # comment signs 
            argPrepend +=  '-s ' + sensors+ ' '
        print('PhoSim added argPrepend is %s'%argPrepend)


        # Do the iteration
        obsId = 9006000
        # opdZkFileName = str("opd.zer" + '.' + testName)
        wfsZkFileName = str("wfs.zer" + '.' + testName)
        opdPssnFileName = "PSSN.txt"
        outputDirName = "pert"
        outputImgDirName = "img"
        iterDefaultDirName = "iter"
        dofInUmFileName = "dofPertInNextIter.mat"
        skyInfoFileName = "skyComCamInfo.txt"
        for iterCount in range(iterNum):

            # Set the observation Id
            phosimCmpt.setSurveyParam(obsId=obsId)

            # The iteration directory
            iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

            # Set the output directory
            outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
            phosimCmpt.setOutputDir(outputDir)
            print('PhoSim outputDir is %s'%outputDir)
            # Set the output image directory
            outputImgDir = os.path.join(baseOutputDir, iterDirName,
                                        outputImgDirName)
            phosimCmpt.setOutputImgDir(outputImgDir)
            print('PhoSim outputImgDir is %s'%outputImgDir)


            # Generate the OPD image
            if genOpd is True:
                
                if selectSensors is 'comcam':
                    argString = phosimCmpt.getComCamOpdArgsAndFilesForPhoSim(
                         cmdSettingFileName=opdCmdSettingsFile)
                    
                elif selectSensors is 'wfs':
                    argString = phosimCmpt.getLsstCamOpdArgsAndFilesForPhoSim(
                        cmdSettingFileName=opdCmdSettingsFile)

                argString = argPrepend + argString
                #argString = '-w $AOCLCOUTPUTPATH ' + argString
                print('Generating OPD with Phosim, argString is \n')
                print(argString)
                phosimCmpt.runPhoSim(argString)

            # Analyze the OPD data
            # Do we need to  analyze the OPD data ? 
            # --> only if need to compare to the WFS results ... 
            # if selectSensors is 'comcam':
            #     phosimCmpt.analyzeComCamOpdData(zkFileName=opdZkFileName,
            #                                 pssnFileName=opdPssnFileName)
            # elif selectSensors is 'wfs':
            #     phosimCmpt.analyzeLsstCamOpdData(zkFileName=opdZkFileName,
            #                                 pssnFileName=opdPssnFileName)

            # Get the PSSN from file
            # pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
            # print("Calculated PSSN is %s." % pssn)

            # # Get the GQ effective FWHM from file
            # gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
            # print("GQ effective FWHM is %.4f." % gqEffFwhm)

            # Set the FWHM data - can't do since they 
            # haven't been calculated ... 
            #listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
            #                                opdPssnFileName, sensorNameList)
            #ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)

            # Use the input sky catalog .... 
            skySim = self._prepareSkySimBySkyFile(inputSkyFilePath)

            # Output the sky information
            skySim, wepCalc = self._outputSkyInfo(outputDir, skyInfoFileName, skySim, wepCalc)

            # Assign the entra- and intra-focal observation Id
            extraObsId = obsId + 1
            intraObsId = obsId + 2

            # Generate the defocal images
            simSeed = 1000
            if selectSensors is 'wfs':
                # I actually don't see here anything specific to comCam.... 
                argStringList = phosimCmpt.getComCamStarArgsAndFilesForPhoSim(
                    extraObsId, intraObsId, skySim, simSeed=simSeed,
                    cmdSettingFileName=comcamCmdSettingsFile,
                instSettingFileName="starSingleExp.inst")

            if genDefocalImg is True:
                for argString in argStringList:
                    #argString = '-w $AOCLCOUTPUTPATH ' + argString
                    argString = argPrepend + argString
                    print('Generating defocal images with Phosim\n')
                    print(argString)
                    phosimCmpt.runPhoSim(argString)

                # Repackage the images based on the image type
                # Again, I don't see here anything comCam - specific...
                # the repackaging calls 
                # /epyc/projects/lsst_comm/phosim_utils/python/lsst/phosim/utils/phosim_repackager.py

                if (isEimg):
                    phosimCmpt.repackageComCamEimgFromPhoSim()
                else:
                    phosimCmpt.repackageComCamAmpImgFromPhoSim()

            # Collect the defocal images
            intraRawExpData = RawExpData()
            intraRawExpDir = os.path.join(outputImgDir,
                                        phosimCmpt.getIntraFocalDirName())
            intraRawExpData.append(intraObsId, 0, intraRawExpDir)

            extraRawExpData = RawExpData()
            extraRawExpDir = os.path.join(outputImgDir,
                                        phosimCmpt.getExtraFocalDirName())
            extraRawExpData.append(extraObsId, 0, extraRawExpDir)

            # Calculate the wavefront error and DOF
            listOfWfErr = wepCalc.calculateWavefrontErrors(
                intraRawExpData, extraRawExpData=extraRawExpData,
                postageImg=postageImg, postageImgDir = postageImgDir)

            # We won't calculate corrections 
            # since we 
            #ofcCalc.calculateCorrections(listOfWfErr)

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
        #pssnFiles = [os.path.join(baseOutputDir, "%s%d" % (iterDefaultDirName, num),
        #            outputImgDirName, opdPssnFileName) for num in range(iterNum)]
        #saveToFilePath = os.path.join(baseOutputDir, "fwhmIters.png")
        #plotFwhmOfIters(pssnFiles, saveToFilePath=saveToFilePath)

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
        chips  = ['00','01','02', 
              '10','11','12',
              '20','21','22']
        rafts = ['00','04', '40', '44']
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
        

    def _preparePhosimCmpt(self, phosimDir, filterType, raInDeg, decInDeg, rotAngInDeg,
                        numPro, isEimg, m1m3ForceError):

        # Set the Telescope facade class
        tele = TeleFacade()
        tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
        tele.setPhoSimDir(phosimDir)

        # Prepare the phosim component
        phosimCmpt = AnalysisPhosimCmpt(tele)

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
                        isEimg,doDeblending, camDimOffset, selectSensors):
        
        if selectSensors is None: # by default
            wepCalc = WEPCalculationFactory.getCalculator(CamType.ComCam, isrDirPath)
        elif selectSensors is 'wfs': # use LsstCam 
            wepCalc = WEPCalculationFactory.getCalculator(CamType.LsstCam, isrDirPath)

        wepCalc.setFilter(filterType)
        wepCalc.setBoresight(raInDeg, decInDeg)
        wepCalc.setRotAng(rotAngInDeg)

        if (isEimg):
            settingFile = wepCalc.getSettingFile()
            settingFile.updateSetting("imageType", "eimage")

        if (doDeblending):
            settingFile = wepCalc.getSettingFile()
            settingFile.updateSetting("doDeblending", "True") 
            
        if camDimOffset  is not None : 
            settingFile = wepCalc.getSettingFile()
            settingFile.updateSetting("camDimOffset", camDimOffset)
            print('camDimOffset is ', settingFile.getSetting("camDimOffset"))


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


