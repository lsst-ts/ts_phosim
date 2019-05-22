import os
import re
import shutil
import numpy as np

from lsst.ts.wep.Utility import runProgram
from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.ctrlIntf.MapSensorNameAndId import MapSensorNameAndId
from lsst.ts.wep.ctrlIntf.SensorWavefrontData import SensorWavefrontData

from lsst.ts.ofc.ctrlIntf.FWHMSensorData import FWHMSensorData

from lsst.ts.phosim.Utility import getConfigDir, sortOpdFileList
from lsst.ts.phosim.OpdMetrology import OpdMetrology


class PhosimCmpt(object):

    def __init__(self, tele):
        """Initialization of PhoSim component class.

        WEP: wavefront estimation pipeline.

        Parameters
        ----------
        tele : TeleFacade
            Telescope instance.
        """

        # Configuration directory
        self.configDir = getConfigDir()

        # Telescope setting file
        settingFilePath = os.path.join(self.configDir,
                                       "phosimCmptSetting.yaml")
        self._phosimCmptSettingFile = ParamReader(filePath=settingFilePath)

        # OPD metrology
        self.metr = OpdMetrology()

        # TeleFacade instance
        self.tele = tele

        # Output directory of data
        self.outputDir = ""

        # Output directory of image
        self.outputImgDir = ""

        # Seed number
        self.seedNum = 0

        # PhoSim parameters
        numPro = self._phosimCmptSettingFile.getSetting("numPro")
        e2ADC = self._phosimCmptSettingFile.getSetting("e2ADC")
        self.phosimParam = {"numPro": numPro,
                            "e2ADC": e2ADC}

    def getNumOfZk(self):
        """Get the number of Zk (annular Zernike polynomial).

        Returns
        -------
        int
            Number of Zk.
        """

        return int(self._phosimCmptSettingFile.getSetting("numOfZk"))

    def getIntraFocalDirName(self):
        """Get the intra-focal directory name.

        Returns
        -------
        str
            Intra-focal directory name.
        """

        return self._phosimCmptSettingFile.getSetting("intraDirName")

    def getExtraFocalDirName(self):
        """Get the extra-focal directory name.

        Returns
        -------
        str
            Extra-focal directory name.
        """

        return self._phosimCmptSettingFile.getSetting("extraDirName")

    def getOpdMetr(self):
        """Get the OPD metrology object.

        OPD: optical path difference.

        Returns
        -------
        OpdMetrology
            OPD metrology object.
        """

        return self.metr

    def setOutputDir(self, outputDir):
        """Set the output directory.

        The output directory will be constructed if there is no existed one.

        Parameters
        ----------
        outputDir : str
            Output directory.
        """

        self._makeDir(outputDir)
        self.outputDir = outputDir

    def _makeDir(self, newDir, exist_ok=True):
        """Make the new directory.

        Super-mkdir; create a leaf directory and all intermediate ones. Works
        like mkdir, except that any intermediate path segment (not just the
        rightmost) will be created if it does not exist.

        Parameters
        ----------
        newDir : str
            New directory.
        exist_ok : bool, optional
            If the target directory already exists, raise an OSError if
            exist_ok is False. Otherwise no exception is raised. (the default
            is True.)
        """

        os.makedirs(newDir, exist_ok=exist_ok)

    def getOutputDir(self):
        """Get the output directory.

        Returns
        -------
        str
            Output directory.
        """

        return self.outputDir

    def setOutputImgDir(self, outputImgDir):
        """Set the output image directory.

        The output image directory will be constructed if there is no existed
        one.

        Parameters
        ----------
        outputImgDir : str
            Output image directory
        """

        self._makeDir(outputImgDir)
        self.outputImgDir = outputImgDir

    def getOutputImgDir(self):
        """Get the output image directory.

        Returns
        -------
        str
            Output image directory
        """

        return self.outputImgDir

    def setSeedNum(self, seedNum):
        """Set the seed number for the M1M3 mirror surface purturbation.

        Parameters
        ----------
        seedNum : int
            Seed number.
        """

        self.seedNum = int(seedNum)

    def getSeedNum(self):
        """Get the seed number for the M1M3 random surface purturbation.

        Returns
        -------
        int or None
            Seed number. None means there is no random purturbation.
        """

        return self.seedNum

    def setPhosimParam(self, numPro, e2ADC):
        """Set the PhoSim simulation parameters.

        Parameters
        ----------
        numPro : int
            Number of processors. The value should be greater than 1.
        e2ADC : int
            Whether to generate amplifier images (1 = true, 0 = false).

        Raises
        ------
        ValueError
            Number of processors should be >= 0.
        ValueError
            e2ADC should be 0 or 1.
        """

        if (numPro > 0):
            self.phosimParam["numPro"] = int(numPro)
        else:
            raise ValueError("Number of processors should be >= 0.")

        if e2ADC in (0, 1):
            self.phosimParam["e2ADC"] = int(e2ADC)
        else:
            raise ValueError("e2ADC should be 0 or 1.")

    def getPhosimParam(self):
        """Get the PhoSim simulation parameters.

        Returns
        -------
        dict
            PhoSim simulation parameters.
        """

        return self.phosimParam

    def setSurveyParam(self, obsId=None, filterType=None, boresight=None,
                       zAngleInDeg=None, rotAngInDeg=None):
        """Set the survey parameters.

        Parameters
        ----------
        obsId : int, optional
            Observation Id. (the default is None.)
        filterType : enum 'FilterType' in lsst.ts.wep.Utility, optional
            Active filter type. (the default is None.)
        boresight : tuple, optional
            Telescope boresight in (ra, decl). (the default is None.)
        zAngleInDeg : float, optional
            Zenith angle in degree. (the default is None.)
        rotAngInDeg : float, optional
            Camera rotation angle in degree between -90 and 90 degrees. (the
            default is None.)
        """

        self.tele.setSurveyParam(obsId=obsId, filterType=filterType,
                                 boresight=boresight, zAngleInDeg=zAngleInDeg,
                                 rotAngInDeg=rotAngInDeg)

    def addOpdFieldXYbyDeg(self, fieldXInDegree, fieldYInDegree):
        """Add the OPD new field X, Y in degree.

        OPD: optical path difference.

        Parameters
        ----------
        fieldXInDegree : float, list, or numpy.ndarray
            New field X in degree.
        fieldYInDegree : float, list, or numpy.ndarray
            New field Y in degree.
        """

        self.metr.addFieldXYbyDeg(fieldXInDegree, fieldYInDegree)

    def accDofInUm(self, dofInUm):
        """Accumulate the aggregated degree of freedom (DOF) in um.

        idx 0-4: M2 dz, dx, dy, rx, ry
        idx 5-9: Cam dz, dx, dy, rx, ry
        idx 10-29: M1M3 20 bending modes
        idx 30-49: M2 20 bending modes

        Parameters
        ----------
        dofInUm : list or numpy.ndarray
            DOF in um.
        """

        self.tele.accDofInUm(dofInUm)

    def setDofInUm(self, dofInUm):
        """Set the accumulated degree of freedom (DOF) in um.

        idx 0-4: M2 dz, dx, dy, rx, ry
        idx 5-9: Cam dz, dx, dy, rx, ry
        idx 10-29: M1M3 20 bending modes
        idx 30-49: M2 20 bending modes

        Parameters
        ----------
        dofInUm : list or numpy.ndarray
            DOF in um.
        """

        self.tele.setDofInUm(dofInUm)

    def getDofInUm(self):
        """Get the accumulated degree of freedom (DOF) in um.

        idx 0-4: M2 dz, dx, dy, rx, ry
        idx 5-9: Cam dz, dx, dy, rx, ry
        idx 10-29: M1M3 20 bending modes
        idx 30-49: M2 20 bending modes

        Returns
        -------
        numpy.ndarray
            DOF in um.
        """

        return self.tele.getDofInUm()

    def saveDofInUmFileForNextIter(self, dofInUm,
                                   dofInUmFileName="dofPertInNextIter.mat"):
        """Save the DOF in um data to file for the next iteration.

        DOF: degree of freedom.

        Parameters
        ----------
        dofInUm : list or numpy.ndarray
            DOF in um.
        dofInUmFileName : str, optional
            File name to save the DOF in um. (the default is
            "dofPertInNextIter.mat".)
        """

        filePath = os.path.join(self.outputDir, dofInUmFileName)
        header = "The followings are the DOF in um:"
        np.savetxt(filePath, np.transpose(dofInUm), header=header)

    def runPhoSim(self, argString):
        """Run the PhoSim program.

        Parameters
        ----------
        argString : str
            Arguments for PhoSim.
        """

        self.tele.runPhoSim(argString)

    def getComCamOpdArgsAndFilesForPhoSim(
            self, cmdFileName="opd.cmd", instFileName="opd.inst",
            logFileName="opdPhoSim.log", cmdSettingFileName="opdDefault.cmd",
            instSettingFileName="opdDefault.inst"):
        """Get the OPD calculation arguments and files of ComCam for the PhoSim
        calculation.

        OPD: optical path difference.
        ComCam: commissioning camera.

        Parameters
        ----------
        cmdFileName : str, optional
            Physical command file name. (the default is "opd.cmd".)
        instFileName : str, optional
            OPD instance file name. (the default is "opd.inst".)
        logFileName : str, optional
            Log file name. (the default is "opdPhoSim.log".)
        cmdSettingFileName : str, optional
            Physical command setting file name. (the default is
            "opdDefault.cmd".)
        instSettingFileName : str, optional
            Instance setting file name. (the default is "opdDefault.inst".)

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        # Set the default ComCam OPD field positions
        self.metr.setDefaultComcamGQ()

        argString = self._getOpdArgsAndFilesForPhoSim(
            cmdFileName, instFileName, logFileName, cmdSettingFileName,
            instSettingFileName)

        return argString

    def _getOpdArgsAndFilesForPhoSim(self, cmdFileName, instFileName,
                                     logFileName, cmdSettingFileName,
                                     instSettingFileName):
        """Get the OPD calculation arguments and files for the PhoSim
        calculation.

        OPD: optical path difference.

        Parameters
        ----------
        cmdFileName : str
            Physical command file name.
        instFileName : str
            OPD instance file name.
        logFileName : str
            Log file name.
        cmdSettingFileName : str
            Physical command setting file name.
        instSettingFileName : str
            Instance setting file name.

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        # Write the command file
        cmdFilePath = self._writePertAndCmdFiles(cmdSettingFileName,
                                                 cmdFileName)

        # Write the instance file
        instSettingFile = self._getInstSettingFilePath(instSettingFileName)
        instFilePath = self.tele.writeOpdInstFile(
            self.outputDir, self.metr, instSettingFile=instSettingFile,
            instFileName=instFileName)

        # Get the argument to run the PhoSim
        argString = self._getPhoSimArgs(logFileName, instFilePath, cmdFilePath)

        return argString

    def _writePertAndCmdFiles(self, cmdSettingFileName, cmdFileName):
        """Write the physical perturbation and command files.

        Parameters
        ----------
        cmdSettingFileName : str
            Physical command setting file name.
        cmdFileName : str
            Physical command file name.

        Returns
        -------
        str
            Command file path.
        """

        # Write the perturbation file
        pertCmdFileName = "pert.cmd"
        pertCmdFilePath = os.path.join(self.outputDir, pertCmdFileName)
        if (not os.path.exists(pertCmdFilePath)):
            self.tele.writePertBaseOnConfigFile(
                self.outputDir, seedNum=self.seedNum, saveResMapFig=True,
                pertCmdFileName=pertCmdFileName)

        # Write the physical command file
        cmdSettingFile = os.path.join(self.configDir, "cmdFile",
                                      cmdSettingFileName)
        cmdFilePath = os.path.join(self.outputDir, cmdFileName)
        if (not os.path.exists(cmdFilePath)):
            self.tele.writeCmdFile(
                self.outputDir, cmdSettingFile=cmdSettingFile,
                pertFilePath=pertCmdFilePath, cmdFileName=cmdFileName)

        return cmdFilePath

    def _getInstSettingFilePath(self, instSettingFileName):
        """Get the instance setting file path.

        Parameters
        ----------
        instSettingFileName : str
            Instance setting file name.

        Returns
        -------
        str
            Instance setting file path.
        """

        instSettingFile = os.path.join(self.configDir, "instFile",
                                       instSettingFileName)

        return instSettingFile

    def _getPhoSimArgs(self, logFileName, instFilePath, cmdFilePath):
        """Get the arguments needed to run the PhoSim.

        Parameters
        ----------
        logFileName : str
            Log file name.
        instFilePath: str
            Instance file path.
        cmdFilePath : str
            Physical command file path.

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        numPro = self.phosimParam["numPro"]
        e2ADC = self.phosimParam["e2ADC"]
        logFilePath = os.path.join(self.outputImgDir, logFileName)

        argString = self.tele.getPhoSimArgs(
            instFilePath, extraCommandFile=cmdFilePath, numPro=numPro,
            outputDir=self.outputImgDir, e2ADC=e2ADC, logFilePath=logFilePath)

        return argString

    def getComCamStarArgsAndFilesForPhoSim(
            self, extraObsId, intraObsId, skySim, simSeed=1000,
            cmdSettingFileName="starDefault.cmd",
            instSettingFileName="starSingleExp.inst"):
        """Get the star calculation arguments and files of ComCam for the
        PhoSim calculation.

        Parameters
        ----------
        extraObsId : int
            Extra-focal observation Id.
        intraObsId : int
            Intra-focal observation Id.
        skySim : SkySim
            Sky simulator
        simSeed : int, optional
            Random number seed. (the default is 1000.)
        cmdSettingFileName : str, optional
            Physical command setting file name. (the default is
            "starDefault.cmd".)
        instSettingFileName : {str}, optional
            Instance setting file name. (the default is "starSingleExp.inst".)

        Returns
        -------
        list[str]
            List of arguments to run the PhoSim.
        """

        # Set the intra- and extra-focal related information
        obsIdList = {"-1": extraObsId,
                     "1": intraObsId}
        instFileNameList = {"-1": "starExtra.inst",
                            "1": "starIntra.inst"}
        logFileNameList = {"-1": "starExtraPhoSim.log",
                           "1": "starIntraPhoSim.log"}

        extraFocalDirName = self.getExtraFocalDirName()
        intraFocalDirName = self.getIntraFocalDirName()
        outImgDirNameList = {"-1": extraFocalDirName,
                             "1": intraFocalDirName}

        # Write the instance and command files of defocal conditions
        cmdFileName = "star.cmd"
        onFocalDofInUm = self.getDofInUm()
        onFocalOutputImgDir = self.outputImgDir
        argStringList = []
        for ii in (-1, 1):

            # Set the observation ID
            self.setSurveyParam(obsId=obsIdList[str(ii)])

            # Camera piston (Change the unit from mm to um)
            pistonInUm = np.zeros(len(onFocalDofInUm))
            pistonInUm[5] = ii * self.tele.getDefocalDistInMm() * 1e3

            # Set the new DOF that considers the piston motion
            self.setDofInUm(onFocalDofInUm + pistonInUm)

            # Update the output image directory
            outputImgDir = os.path.join(onFocalOutputImgDir,
                                        outImgDirNameList[str(ii)])
            self.setOutputImgDir(outputImgDir)

            # Get the argument to run the phosim
            argString = self.getStarArgsAndFilesForPhoSim(
                skySim, cmdFileName=cmdFileName,
                instFileName=instFileNameList[str(ii)],
                logFileName=logFileNameList[str(ii)], simSeed=simSeed,
                cmdSettingFileName=cmdSettingFileName,
                instSettingFileName=instSettingFileName)
            argStringList.append(argString)

        # Put the internal state back to the focal plane condition
        self.setDofInUm(onFocalDofInUm)
        self.setOutputImgDir(onFocalOutputImgDir)

        return argStringList

    def getStarArgsAndFilesForPhoSim(self, skySim,
                                     cmdFileName="star.cmd",
                                     instFileName="star.inst",
                                     logFileName="starPhoSim.log",
                                     simSeed=1000,
                                     cmdSettingFileName="starDefault.cmd",
                                     instSettingFileName="starSingleExp.inst"):
        """Get the star calculation arguments and files for the PhoSim
        calculation.

        Parameters
        ----------
        skySim : SkySim
            Sky simulator
        cmdFileName : str, optional
            Physical command file name. (the default is "star.cmd".)
        instFileName : str, optional
            Star instance file name. (the default is "star.inst".)
        logFileName : str, optional
            Log file name. (the default is "starPhoSim.log".)
        simSeed : int, optional
            Random number seed. (the default is 1000)
        cmdSettingFileName : str, optional
            Physical command setting file name. (the default is
            "starDefault.cmd".)
        instSettingFileName : str, optional
            Instance setting file name. (the default is "starSingleExp.inst".)

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        # Write the command file
        cmdFilePath = self._writePertAndCmdFiles(cmdSettingFileName,
                                                 cmdFileName)

        # Write the instance file
        instSettingFile = self._getInstSettingFilePath(instSettingFileName)
        instFilePath = self.tele.writeStarInstFile(
            self.outputDir, skySim, simSeed=simSeed, sedName="sed_flat.txt",
            instSettingFile=instSettingFile, instFileName=instFileName)

        # Get the argument to run the PhoSim
        argString = self._getPhoSimArgs(logFileName, instFilePath, cmdFilePath)

        return argString

    def analyzeComCamOpdData(self, zkFileName="opd.zer",
                             pssnFileName="PSSN.txt"):
        """Analyze the ComCam OPD data.

        ComCam: Commissioning camera.
        OPD: optical path difference.
        PSSN: normalized point source sensitivity.

        Parameters
        ----------
        zkFileName : str, optional
            OPD in zk file name. (the default is "opd.zer".)
        pssnFileName : str, optional
            PSSN file name. (the default is "PSSN.txt".)
        """

        self._writeOpdZkFile(zkFileName)
        self._writeOpdPssnFile(pssnFileName)

    def _writeOpdZkFile(self, zkFileName):
        """Write the OPD in zk file.

        OPD: optical path difference.

        Parameters
        ----------
        zkFileName : str
            OPD in zk file name.
        """

        filePath = os.path.join(self.outputImgDir, zkFileName)
        opdData = self._mapOpdToZk()
        header = "The followings are OPD in um from z4 to z22:"
        np.savetxt(filePath, opdData, header=header)

    def _mapOpdToZk(self):
        """Map the OPD to the basis of annular Zernike polynomial (Zk).

        OPD: optical path difference.

        Returns
        -------
        numpy.ndarray
            Zk data from OPD. This is a 2D array. The row is the OPD index and
            the column is z4 to z22 in um. The order of OPD index is based on
            the file name.
        """

        # Get the sorted OPD file list
        opdFileList = self._getOpdFileInDir(self.outputImgDir)

        # Map the OPD to the Zk basis and do the collection
        numOfZk = self.getNumOfZk()
        opdData = np.zeros((len(opdFileList), numOfZk))
        for idx, opdFile in enumerate(opdFileList):

            # z1 to z22 (22 terms)
            zk = self.metr.getZkFromOpd(opdFitsFile=opdFile)[0]

            # Only need to collect z4 to z22
            initIdx = 3
            opdData[idx, :] = zk[initIdx:initIdx + numOfZk]

        return opdData

    def _getOpdFileInDir(self, opdDir):
        """Get the sorted OPD files in the directory.

        OPD: Optical path difference.

        Parameters
        ----------
        opdDir : str
            OPD file directory.

        Returns
        -------
        list
            List of sorted OPD files.
        """

        # Get the files
        opdFileList = []
        fileList = self._getFileInDir(opdDir)
        for file in fileList:
            fileName = os.path.basename(file)
            m = re.match(r"\Aopd_\d+_(\d+).fits.gz", fileName)
            if (m is not None):
                opdFileList.append(file)

        # Do the sorting of file name
        sortedOpdFileList = sortOpdFileList(opdFileList)

        return sortedOpdFileList

    def _getFileInDir(self, fileDir):
        """Get the files in the directory.

        Parameters
        ----------
        fileDir : str
            File directory.

        Returns
        -------
        list
            List of files.
        """

        fileList = []
        for name in os.listdir(fileDir):
            filePath = os.path.join(fileDir, name)
            if os.path.isfile(filePath):
                fileList.append(filePath)

        return fileList

    def _writeOpdPssnFile(self, pssnFileName):
        """Write the OPD PSSN in file.

        OPD: optical path difference.
        PSSN: normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.
        """

        filePath = os.path.join(self.outputImgDir, pssnFileName)

        # Calculate the PSSN
        pssnList, gqEffPssn = self._calcComCamOpdPssn()

        # Calculate the FWHM
        effFwhmList, gqEffFwhm = self._calcComCamOpdEffFwhm(pssnList)

        # Append the list to write the data into file
        pssnList.append(gqEffPssn)
        effFwhmList.append(gqEffFwhm)

        # Stack the data
        data = np.vstack((pssnList, effFwhmList))

        # Write to file
        header = "The followings are PSSN and FWHM (in arcsec) data. The final number is the GQ value."
        np.savetxt(filePath, data, header=header)

    def _calcComCamOpdPssn(self):
        """Calculate the ComCam PSSN of OPD.

        ComCam: commissioning camera.
        OPD: optical path difference.
        PSSN: normalized point source sensitivity.
        GQ: gaussian quadrature.

        Returns
        -------
        list
            PSSN list.
        float
            GQ effective PSSN.
        """

        opdFileList = self._getOpdFileInDir(self.outputImgDir)

        wavelengthInUm = self.tele.getRefWaveLength() * 1e-3
        pssnList = []
        for opdFile in opdFileList:
            pssn = self.metr.calcPSSN(wavelengthInUm, opdFitsFile=opdFile)
            pssnList.append(pssn)

        # Calculate the GQ effectice PSSN
        self._setComCamWgtRatio()
        gqEffPssn = self.metr.calcGQvalue(pssnList)

        return pssnList, gqEffPssn

    def _setComCamWgtRatio(self):
        """Set the ComCam weighting ratio.

        ComCam: Commissioning camera.
        """

        comcamWtRatio = np.ones(9)
        self.metr.setWeightingRatio(comcamWtRatio)

    def _calcComCamOpdEffFwhm(self, pssnList):
        """Calculate the ComCam effective FWHM of OPD.

        ComCam: commissioning camera.
        FWHM: full width and half maximum.
        PSSN: normalized point source sensitivity.
        GQ: Gaussian quadrature.

        Parameters
        ----------
        pssnList : list
            List of PSSN.

        Returns
        -------
        list
            Effective FWHM list.
        float
            GQ effective FWHM of ComCam.
        """

        # Calculate the list of effective FWHM
        effFwhmList = []
        for pssn in pssnList:
            effFwhm = self.metr.calcFWHMeff(pssn)
            effFwhmList.append(effFwhm)

        # Calculate the GQ effectice FWHM
        self._setComCamWgtRatio()
        gqEffFwhm = self.metr.calcGQvalue(effFwhmList)

        return effFwhmList, gqEffFwhm

    def mapOpdDataToListOfWfErr(self, opdZkFileName, refSensorNameList):
        """Map the OPD data to the list of wavefront error.

        OPD: Optical path difference.

        Parameters
        ----------
        opdZkFileName : str
            OPD zk file name.
        refSensorNameList : list
            Reference sensor name list.

        Returns
        -------
        list[SensorWavefrontData]
            List of SensorWavefrontData object.
        """

        opdZk = self._getZkFromFile(opdZkFileName)

        mapSensorNameAndId = MapSensorNameAndId()
        sensorIdList = mapSensorNameAndId.mapSensorNameToId(refSensorNameList)

        listOfWfErr = []
        for sensorId, zk in zip(sensorIdList, opdZk):

            sensorWavefrontData = SensorWavefrontData()
            sensorWavefrontData.setSensorId(sensorId)
            sensorWavefrontData.setAnnularZernikePoly(zk)

            listOfWfErr.append(sensorWavefrontData)

        return listOfWfErr

    def _getZkFromFile(self, zkFileName):
        """Get the zk (z4-z22) from file.

        Parameters
        ----------
        zkFileName : str
            Zk file name.

        Returns
        -------
        numpy.ndarray
            zk matrix. The colunm is z4-z22. The raw is each data point.
        """

        filePath = os.path.join(self.outputImgDir, zkFileName)
        zk = np.loadtxt(filePath)

        return zk

    def getOpdPssnFromFile(self, pssnFileName):
        """Get the OPD PSSN from file.

        OPD: optical path difference.
        PSSN: normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.

        Returns
        -------
        numpy.ndarray
            PSSN.
        """

        data = self._getDataOfPssnFile(pssnFileName)
        pssn = data[0, :-1]

        return pssn

    def _getDataOfPssnFile(self, pssnFileName):
        """Get the data of the PSSN file.

        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.

        Returns
        -------
        numpy.ndarray
            Data of the PSSN file.
        """

        filePath = os.path.join(self.outputImgDir, pssnFileName)
        data = np.loadtxt(filePath)

        return data

    def getOpdGqEffFwhmFromFile(self, pssnFileName):
        """Get the OPD GQ effective FWHM from file.

        OPD: Optical path difference.
        GQ: Gaussian quadrature.
        FWHM: Full width at half maximum.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.

        Returns
        -------
        float
            OPD GQ effective FWHM.
        """

        data = self._getDataOfPssnFile(pssnFileName)
        gqEffFwhm = data[1, -1]

        return gqEffFwhm

    def getListOfFwhmSensorData(self, pssnFileName, refSensorNameList):
        """Get the list of FWHM sensor data based on the OPD PSSN file.

        FWHM: Full width at half maximum.
        OPD: Optical path difference.
        PSSN: Normalized point source sensitivity.

        Parameters
        ----------
        pssnFileName : str
            PSSN file name.
        refSensorNameList : list
            Reference sensor name list.

        Returns
        -------
        list[FWHMSensorData]
            List of FWHMSensorData which contains the sensor Id and FWHM data.
        """

        # Get the FWHM data from the PSSN file
        # The first row is the PSSN and the second one is the FWHM
        # The final element in each row is the GQ value
        data = self._getDataOfPssnFile(pssnFileName)
        fwhmData = data[1, :-1]

        mapSensorNameAndId = MapSensorNameAndId()
        sensorIdList = mapSensorNameAndId.mapSensorNameToId(refSensorNameList)

        listOfFWHMSensorData = []
        for sensorId, fwhm in zip(sensorIdList, fwhmData):
            fwhmSensorData = FWHMSensorData(sensorId, np.array([fwhm]))
            listOfFWHMSensorData.append(fwhmSensorData)

        return listOfFWHMSensorData

    def repackageComCamImgFromPhoSim(self):
        """Repackage the ComCam amplifier images from PhoSim to the single 16
        extension MEFs for processing.

        ComCam: commissioning camera.
        MEF: multi-extension frames.
        """

        # Make a temp directory
        tmpDirPath = os.path.join(self.outputImgDir, "tmp")
        self._makeDir(tmpDirPath)

        intraFocalDirName = self.getIntraFocalDirName()
        extraFocalDirName = self.getExtraFocalDirName()
        for imgType in (intraFocalDirName, extraFocalDirName):

            # Repackage the images to that temp directory
            command = "phosim_repackager.py"
            phosimAmgImgDir = os.path.join(self.outputImgDir, imgType)
            argstring = "%s --out_dir=%s" % (phosimAmgImgDir, tmpDirPath)
            runProgram(command, argstring=argstring)

            # Remove the image data in the original directory
            argString = "-rf %s/lsst_a_*.fits.gz" % phosimAmgImgDir
            runProgram("rm", argstring=argString)

            # Put the repackaged data into the image directory
            argstring = "%s/*.fits %s" % (tmpDirPath, phosimAmgImgDir)
            runProgram("mv", argstring=argstring)

        # Remove the temp directory
        shutil.rmtree(tmpDirPath)

    def reorderAndSaveWfErrFile(self, listOfWfErr, refSensorNameList,
                                zkFileName="wfs.zer"):
        """Reorder the wavefront error in the wavefront error list according to
        the reference sensor name list and save to a file.

        The unexisted wavefront error will be a numpy zero array. The unit is
        um.

        Parameters
        ----------
        listOfWfErr : list[SensorWavefrontData]
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
        filePath = os.path.join(self.outputImgDir, zkFileName)
        wfsData = self._getWfErrValuesAndStackToMatrix(reorderedWfErrMap)
        header = "The followings are ZK in um from z4 to z22:"
        np.savetxt(filePath, wfsData, header=header)

    def _transListOfWfErrToMap(self, listOfWfErr):
        """Transform the list of wavefront error to map.

        Parameters
        ----------
        listOfWfErr : list[SensorWavefrontData]
            List of SensorWavefrontData object.

        Returns
        -------
        dict
            Calculated wavefront error. The dictionary key [str] is the
            abbreviated sensor name (e.g. R22_S11). The dictionary item
            [numpy.ndarray] is the averaged wavefront error (z4-z22) in um.
        """

        mapSensorNameAndId = MapSensorNameAndId()

        wfErrMap = dict()
        for sensorWavefrontData in listOfWfErr:
            sensorId = sensorWavefrontData.getSensorId()
            sensorNameList = mapSensorNameAndId.mapSensorIdToName(sensorId)[0]
            sensorName = sensorNameList[0]

            avgErrInUm = sensorWavefrontData.getAnnularZernikePoly()

            wfErrMap[sensorName] = avgErrInUm

        return wfErrMap

    def _getWfErrValuesAndStackToMatrix(self, wfErrMap):
        """Get the wavefront errors and stack them to be a matrix.

        Parameters
        ----------
        wfErrMap : dict
            Calculated wavefront error. The dictionary key [str] is the
            abbreviated sensor name (e.g. R22_S11). The dictionary item
            [numpy.ndarray] is the averaged wavefront error (z4-z22) in um.

        Returns
        -------
        numpy.ndarray
            Wavefront errors as a matrix. The column is z4-z22 in um. The row
            is the individual sensor. The order is the same as the input of
            wfErrMap.
        """

        numOfZk = self.getNumOfZk()
        valueMatrix = np.empty((0, numOfZk))
        for wfErr in wfErrMap.values():
            valueMatrix = np.vstack((valueMatrix, wfErr))

        return valueMatrix


if __name__ == "__main__":
    pass
