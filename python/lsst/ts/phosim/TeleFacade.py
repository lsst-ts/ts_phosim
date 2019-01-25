import os
import re
import numpy as np

from lsst.sims.utils import ObservationMetaData
from lsst.ts.wep.Utility import FilterType, mapFilterRefToG

from lsst.ts.phosim.CamSim import CamSim
from lsst.ts.phosim.M1M3Sim import M1M3Sim
from lsst.ts.phosim.M2Sim import M2Sim
from lsst.ts.phosim.PhosimCommu import PhosimCommu
from lsst.ts.phosim.Utility import SurfaceType, CamDistType, mapSurfNameToEnum, \
    createObservation


class TeleFacade(object):

    # Reference wavelength of reference filter (FilterType.REF)
    WAVELENGTH_IN_NM = 500

    def __init__(self, configFilePath=None):
        """Initialization of telescope facade class.

        This class uses the facade pattern that the high level class telescope
        helps to write the perturbations of camera, M1M3, and M2 into the
        PhoSim by the interface to PhoSim.

        Parameters
        ----------
        configFilePath : str, optional
            Telescope configuration file path. (the default is None.)
        """

        self.cam = None
        self.M1M3 = None
        self.M2 = None
        self.phoSimCommu = PhosimCommu()

        self.dofInUm = np.zeros(50)

        self.instName = "lsst"
        self.defocalDisInMm = 1.5

        # Default Observation
        self.obs = createObservation()

        self.sensorOn = {"sciSensorOn": True,
                         "wfSensorOn": True,
                         "guidSensorOn": False}

        self.configFile = configFilePath

    def getDefocalDisInMm(self):
        """Get the defocal distance in mm.

        Returns
        -------
        float
            defocal distance in mm.
        """

        return self.defocalDisInMm

    def setObservation(self, obs):
        """Set the survey parameters from an observation.

        Parameters
        ----------
        obs : ObservationMetaData
            The observation parameters.
        """
        self.obs = obs

    def getObservation(self):
        """Get the observation.

        Returns
        -------
        ObservationMetaData
            The observation parameters.
        """
        return self.obs

    def setSensorOn(self, sciSensorOn=True, wfSensorOn=True,
                    guidSensorOn=False):
        """Set the sensor on.

        Parameters
        ----------
        sciSensorOn : bool, optional
            Scientific sensors are on. (the default is True.)
        wfSensorOn : bool, optional
            Wavefront sensors are on. (the default is True.)
        guidSensorOn : bool, optional
            Guider sensors are on. (the default is False.)
        """

        if (self._isBool(sciSensorOn)):
            self.sensorOn["sciSensorOn"] = sciSensorOn

        if (self._isBool(wfSensorOn)):
            self.sensorOn["wfSensorOn"] = wfSensorOn

        if (self._isBool(guidSensorOn)):
            self.sensorOn["guidSensorOn"] = guidSensorOn

    def _isBool(self, input):
        """Check the input is boolean type or not.

        Parameters
        ----------
        input : bool
            Input.

        Returns
        -------
        bool
            True if the type of input is boolean.
        """

        if isinstance(input, bool):
            return True
        else:
            return False

    def setConfigFile(self, configFilePath):
        """Set the telescope configuration file.

        Parameters
        ----------
        configFilePath : str
            Configuration file path.
        """

        self.configFile = configFilePath

    def getConfigValue(self, varName, index=1):
        """Get the value of certain variable defined in the configuration file.

        Parameters
        ----------
        varName : str
            Name of variable.
        index : int, optional
            Index of value. (the default is 1.)

        Returns
        -------
        float, int, or str
            Variable value.
        """

        # Read the file
        fid = open(self.configFile)

        # Search for the value of certain variable
        value = None
        for line in fid:

            # Strip the line
            line = line.strip()

            # Get the element of line
            lineArray = line.split()

            # Get the value
            if (len(lineArray) > 0) and (lineArray[0] == varName):
                value = lineArray[index]
                break

        # Close the file
        fid.close()

        # Change the value type if necessary
        try:
            value = float(value)
            if (value == int(value)):
                value = int(value)
        except Exception as ValueError:
            pass

        return value

    def runPhoSim(self, argString):
        """

        Run the PhoSim program.

        Arguments:
            argString {[str]} -- Arguments for PhoSim.
        """

        self.phoSimCommu.runPhoSim(argstring=argString)

    def getPhoSimArgs(self, instFilePath, extraCommandFile=None, numPro=1,
                      numThread=1, outputDir=None, sensorName=None,
                      e2ADC=1, logFilePath=None):
        """Get the arguments needed to run the PhoSim.

        Parameters
        ----------
        instanceFile : str
            Instance catalog file.
        extraCommandFile : str, optional
            Command file to modify the default physics. (the default is None.)
        numProc : int, optional
             Number of processors. (the default is 1.)
        numThread : int, optional
            Number of threads. (the default is 1.)
        outputDir : str, optional
            Output image directory. (the default is None.)
        sensorName : str, optional
            Sensor chip specification (e.g., all, R22_S11, "R22_S11|R22_S12").
            (the default is None.)
        e2ADC : int, optional
            Whether to generate amplifier images (1 = true, 0 = false). (the
            default is 1.)
        logFilePath : str, optional
            Log file path of PhoSim calculation. (the default is None.)

        Returns
        -------
        str
            Arguments to run the PhoSim.
        """

        argString = self.phoSimCommu.getPhoSimArgs(
            instFilePath, extraCommandFile=extraCommandFile, numProc=numPro,
            numThread=numThread, outputDir=outputDir, instrument=self.instName,
            sensorName=sensorName, e2ADC=e2ADC, logFilePath=logFilePath)

        return argString

    def setInstName(self, instName):
        """ Set the instrument name from the instrument file.

        Parameters
        ----------
        instName : str
            Instrument name (e.g. "lsst", "comcam10").

        Raises
        ------
        ValueError
            Cannot get the instrument name.
        """

        # Get the instrument name
        m = re.match(r"([a-z]+)(?:(\d+))?$", instName)
        if (m is None):
            raise ValueError("Cannot get the instrument name: %s." % instName)
        instName = m.groups()[0]

        # Decide the defocal distance offset in mm
        defocalOffset = m.groups()[1]
        if (defocalOffset is not None):
            defocalOffset = float(defocalOffset) / 10
        else:
            # Default defocal distance is 1.5 mm
            defocalOffset = 1.5

        self.instName = instName
        self.defocalDisInMm = defocalOffset

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

        self.dofInUm = np.array(dofInUm, dtype=float)

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

        self.dofInUm += np.array(dofInUm, dtype=float)

    def setSubSysConfigDir(self, camDataDir=None, M1M3dataDir=None,
                           M2dataDir=None, phosimDir=None):
        """Set the subsystem data directory.

        Parameters
        ----------
        camDataDir : str, optional
            Camera data directory. (the default is None.)
        M1M3dataDir : str, optional
            M1M3 data directory. (the default is None.)
        M2dataDir : str, optional
            M2 data directory. (the default is None.)
        phosimDir : str, optional
            PhoSim directory. (the default is None.)
        """

        # Get the number of zernike terms used in Zemax
        numTerms = self.getConfigValue("znPert")

        if (camDataDir is not None):
            self.cam = CamSim()
            self.cam.setCamDataDir(camDataDir)

        if (M1M3dataDir is not None):
            self.M1M3 = M1M3Sim()
            self.M1M3.setMirrorDataDir(M1M3dataDir)
            self.M1M3.config(numTerms=numTerms)

        if (M2dataDir is not None):
            self.M2 = M2Sim()
            self.M2.setMirrorDataDir(M2dataDir)
            self.M2.config(numTerms=numTerms)

        if (phosimDir is not None):
            self.phoSimCommu.setPhoSimDir(phosimDir)

    def writeAccDofFile(self, outputFileDir, dofFileName="pert.mat"):
        """Write the accumulated degree of freedom (DOF) in um to file.

        Parameters
        ----------
        outputFileDir : str
            Output file directory.
        dofFileName : str, optional
            DOF file name. (the default is "pert.mat".)

        Returns
        -------
        str
            DOF file path.
        """

        dofFilePath = os.path.join(outputFileDir, dofFileName)
        np.savetxt(dofFilePath, self.dofInUm)

        return dofFilePath

    def writeCmdFile(self, cmdFileDir, cmdSettingFile=None, pertFilePath=None,
                     cmdFileName="taskPert.cmd"):
        """Write the physical command file.

        OPD: Optical path difference.

        Parameters
        ----------
        cmdFileDir : str
            Directory to the OPD command file.
        cmdSettingFile : str, optional
            Physical command setting file path. (the default is None.)
        pertFilePath : str, optional
            Subsystem perturbation command file path. (the default is None.)
        cmdFileName : str, optional
            Command file name. (the default is "taskPert.cmd".)

        Returns
        -------
        str
            Command file path.
        """

        # Command file path
        cmdFilePath = os.path.join(cmdFileDir, cmdFileName)
        self.phoSimCommu.writeToFile(cmdFilePath, content="", mode="w")

        # Write the physical setting
        if (cmdSettingFile is not None):
            self.phoSimCommu.writeToFile(cmdFilePath,
                                         sourceFile=cmdSettingFile)

        # Add the subsystem perturbation
        if (pertFilePath is not None):
            self.phoSimCommu.writeToFile(cmdFilePath, sourceFile=pertFilePath)

        return cmdFilePath

    def writeStarInstFile(self, instFileDir, skySim, simSeed=1000,
                          sedName="sed_flat.txt", instOverrideFile=None,
                          instFileName="star.inst"):
        """Write the star instance file.

        Parameters
        ----------
        instFileDir : str
            Directory to instance file.
        skySim : SkySim
            SkySim object.
        simSeed : int, optional
            Random number seed. (the default is 1000.)
        sedName : str, optional
            The name of the SED file with a file path that is relative to the
            data directory in PhoSim. (the default is "sed_flat.txt".)
        instOverrideFile : str optional
            Instance setting file to with PhoSim arguments to override observation.
            (the default is None.)
        instFileName : str, optional
            Star instance file name. (the default is "star.inst".)

        Returns
        -------
        str
            Instance file path.
        """

        # Instance file path
        instFilePath = os.path.join(instFileDir, instFileName)

        # Write the default instance setting
        obsId = self.obs.OpsimMetaData["obsHistID"]

        ra = self.obs.pointingRA
        dec = self.obs.pointingDec

        rot = self.obs.rotSkyPos
        mjd = self.obs.mjd.TAI

        # Handle with FilterType.REF
        filt = mapFilterRefToG(FilterType.fromString(self.obs.bandpass))
        self.obs.setBandpassM5andSeeing(filt.toString())

        self.phoSimCommu.writeObsHeader(instFilePath, self.obs, mode="w")

        if (instFilePath is not None):
            self.phoSimCommu.writeToFile(instFilePath,
                                         sourceFile=instOverrideFile)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Set the camera configuration
        sciSensorOn = self.sensorOn["sciSensorOn"]
        wfSensorOn = self.sensorOn["wfSensorOn"]
        guidSensorOn = self.sensorOn["guidSensorOn"]
        content += self.phoSimCommu.doCameraConfig(
            sciSensorOn=sciSensorOn, wfSensorOn=wfSensorOn,
            guidSensorOn=guidSensorOn)

        # Use the SED file of single wavelendth if the reference filter is used.
        filterType = FilterType.fromString(self.obs.bandpass)
        if (filterType == FilterType.REF):
            self._writeSedFileIfPhoSimDirSet()
            sedName = "sed_%s.txt" % int(self.WAVELENGTH_IN_NM)

        # Write the star source
        for ii in range(len(skySim.starId)):
            content += self.phoSimCommu.generateStar(
                skySim.starId[ii], skySim.ra[ii], skySim.decl[ii],
                skySim.mag[ii], sedName)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        return instFilePath

    def _getFilterIdInPhoSim(self):
        """Get the active filter Id used in PhoSim.

        Returns
        -------
        int
            Active filter ID in PhoSim.
        """

        filterType = FilterType.fromString(self.obs.bandpass)
        mappedFilterType = mapFilterRefToG(filterType)

        return self.phoSimCommu.getFilterId(mappedFilterType)

    def writeOpdInstFile(self, instFileDir, opdMetr, instSettingFile=None,
                         instFileName="opd.inst"):
        """Write the optical path difference (OPD) instance file.

        Parameters
        ----------
        instFileDir : str
            Directory to instance file.
        opdMetr : OpdMetrology
            OpdMetrology object.
        instSettingFile : str, optional
            Instance setting file. (the default is None.)
        instFileName : str, optional
            OPD instance file name. (the default is "opd.inst".)

        Returns
        -------
        str
            Instance file path.
        """

        # Instance file path
        instFilePath = os.path.join(instFileDir, instFileName)

        # Get the observation ID
        obsId = self.obs.OpsimMetaData["obsHistID"]

        # Get the filter ID in PhoSim
        aFilterId = self._getFilterIdInPhoSim()

        # Add the sky information
        ra = self.obs.pointingRA
        dec = self.obs.pointingDec

        rot = self.obs.rotSkyPos

        # Write the default instance setting
        self.phoSimCommu.getOpdInstance(obsId, aFilterId, ra=ra, dec=dec, 
                                        rot=rot, filePath=instFilePath)
        if (instSettingFile is not None):
            self.phoSimCommu.writeToFile(instFilePath,
                                         sourceFile=instSettingFile)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Write the OPD source
        for ii in range(len(opdMetr.fieldX)):
            content += self.phoSimCommu.generateOpd(
                ii, opdMetr.fieldX[ii], opdMetr.fieldY[ii],
                self.WAVELENGTH_IN_NM)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        # Write the OPD SED file if necessary
        self._writeSedFileIfPhoSimDirSet()

        return instFilePath

    def _writeSedFileIfPhoSimDirSet(self):
        """Write the SED file if the PhoSim directory is set.

        SED: Spectral energy distribution.
        """

        if os.path.isdir(self.phoSimCommu.phosimDir):
            self.phoSimCommu.writeSedFile(self.WAVELENGTH_IN_NM)
        else:
            print("Do not inspect the SED file for no PhoSim directory.")

    def writePertBaseOnConfigFile(self, pertCmdFileDir, seedNum=None,
                                  M1M3ForceError=0.05, saveResMapFig=False,
                                  pertCmdFileName="pert.cmd"):
        """Write the perturbation command file based on the telescope
        configuration file.

        Parameters
        ----------
        pertCmdFileDir : str
            Directory to the pertubation command file.
        seedNum : int, optional
            Random seed number. If the value is not None, the M1M3 mirror
            will generate a random surface error. (the default is None.)
        M1M3ForceError : float, optional
            Ratio of actuator force error. (the default is 0.05.)
        saveResMapFig : bool, optional
            Save the mirror surface residue map or not. (the default is False.)
        pertCmdFileName : str, optional
            Perturbation command file name. (the default is "pert.cmd".)

        Returns
        -------
        str
            Perturbation command file path.
        """

        # Mirror surface residue file name
        M1resFileName = "M1res.txt"
        M2resFileName = "M2res.txt"
        M3resFileName = "M3res.txt"

        # Mirror Zc file
        M1M3zcFileName = "M1M3zlist.txt"
        M2zcFileName = "M2zlist.txt"

        # Path of perturbation command file
        pertCmdFilePath = os.path.join(pertCmdFileDir, pertCmdFileName)
        M1resFilePath = os.path.join(pertCmdFileDir, M1resFileName)
        M2resFilePath = os.path.join(pertCmdFileDir, M2resFileName)
        M3resFilePath = os.path.join(pertCmdFileDir, M3resFileName)
        M1M3zcFilePath = os.path.join(pertCmdFileDir, M1M3zcFileName)
        M2zcFilePath = os.path.join(pertCmdFileDir, M2zcFileName)

        # Perturbation command
        content = ""

        # Get the zenith angle in radian
        zAngleInDeg = 90 - self.obs.OpsimMetaData["altitude"]
        zAngleInRad = np.deg2rad(zAngleInDeg)

        # Get the numeber of grid used in Zemax
        surfaceGridN = self.getConfigValue("surfaceGridN")

        # Write the camera perturbation command file
        if (self.M1M3 is not None):

            # Do the gravity correction
            printthzInM = self.M1M3.getPrintthz(zAngleInRad)

            # Add the surface error if necessary
            randSurfInM = None
            if (seedNum is not None):
                randSurfInM = self.M1M3.genMirSurfRandErr(
                    zAngleInRad, M1M3ForceError=M1M3ForceError,
                    seedNum=seedNum)

            # Do the temperature correction
            M1M3TBulk = self.getConfigValue("M1M3TBulk")
            M1M3TxGrad = self.getConfigValue("M1M3TxGrad")
            M1M3TyGrad = self.getConfigValue("M1M3TyGrad")
            M1M3TzGrad = self.getConfigValue("M1M3TzGrad")
            M1M3TrGrad = self.getConfigValue("M1M3TrGrad")
            tempCorrInUm = self.M1M3.getTempCorr(M1M3TBulk, M1M3TxGrad,
                                                 M1M3TyGrad, M1M3TzGrad,
                                                 M1M3TrGrad)

            # Set the mirror surface in mm
            if (randSurfInM is not None):
                mirrorSurfInUm = (printthzInM + randSurfInM) * 1e6 + \
                    tempCorrInUm
            else:
                mirrorSurfInUm = printthzInM * 1e6 + tempCorrInUm
            self.M1M3.setSurfAlongZ(mirrorSurfInUm)

            resFile = [M1resFilePath, M3resFilePath]
            self.M1M3.writeMirZkAndGridResInZemax(
                resFile=resFile, surfaceGridN=surfaceGridN,
                writeZcInMnToFilePath=M1M3zcFilePath)

            # Get the Zk in mm
            zkInMm = np.loadtxt(M1M3zcFilePath)

            # Do the surface perturbation
            surfList = [SurfaceType.M1, SurfaceType.M3]
            surfIdList = []
            for ii in range(len(surfList)):
                surf = surfList[ii]
                surfId = self.phoSimCommu.getSurfaceId(surf)
                content += self.phoSimCommu.doSurfPert(surfId, zkInMm)

                # Collect the surface ID
                surfIdList.append(surfId)

                # Do the surface residue map perturbation
                content += self.phoSimCommu.doSurfMapPert(surfId,
                                                          resFile[ii], 1)

            # Do the surface linkage
            content += self.phoSimCommu.doSurfLink(surfIdList[1],
                                                   surfIdList[0])

        if (self.M2 is not None):

            # Do the gravity correction
            printthzInUm = self.M2.getPrintthz(zAngleInRad)

            # Do the temperature correction
            M2TzGrad = self.getConfigValue("M2TzGrad")
            M2TrGrad = self.getConfigValue("M2TrGrad")
            tempCorrInUm = self.M2.getTempCorr(M2TzGrad, M2TrGrad)

            # Set the mirror surface in mm
            mirrorSurfInUm = printthzInUm + tempCorrInUm
            self.M2.setSurfAlongZ(mirrorSurfInUm)
            self.M2.writeMirZkAndGridResInZemax(
                resFile=M2resFilePath, surfaceGridN=surfaceGridN,
                writeZcInMnToFilePath=M2zcFilePath)

            # Get the Zk in mm
            zkInMm = np.loadtxt(M2zcFilePath)

            # Do the surface perturbation
            surfId = self.phoSimCommu.getSurfaceId(SurfaceType.M2)
            content += self.phoSimCommu.doSurfPert(surfId, zkInMm)

            # Do the surface residue map perturbation
            content += self.phoSimCommu.doSurfMapPert(surfId, M2resFilePath, 1)

        if (self.cam is not None):
            # Set the camera rotation angle
            rotAngInDeg = self.obs.rotSkyPos
            self.cam.setRotAngInDeg(rotAngInDeg)

            # Set the temperature information
            tempInDegC = self.getConfigValue("camTB")
            self.cam.setBodyTempInDegC(tempInDegC)

            # Write the perturbation file
            for distType in CamDistType:
                # Get the surface ID
                surfaceType = self._getPhoSimCamSurf(distType.name)
                surfId = self.phoSimCommu.getSurfaceId(surfaceType)

                # Do the perturbation
                zkInMm = self.cam.getCamDistortionInMm(zAngleInRad, distType)
                content += self.phoSimCommu.doSurfPert(surfId, zkInMm)

        # Write the perturbation command to file
        self.phoSimCommu.writeToFile(pertCmdFilePath, content=content,
                                     mode="w")

        # Save the mirror residue map if necessary
        if (saveResMapFig):
            if (self.M1M3 is not None):
                resFile = [M1resFilePath, M3resFilePath]
                writeToResMapFilePath1 = os.path.splitext(M1resFilePath)[0] + \
                    ".png"
                writeToResMapFilePath3 = os.path.splitext(M3resFilePath)[0] + \
                    ".png"
                writeToResMapFilePath = [writeToResMapFilePath1,
                                         writeToResMapFilePath3]
                self.M1M3.showMirResMap(
                    resFile=resFile,
                    writeToResMapFilePath=writeToResMapFilePath)

            if (self.M2 is not None):
                writeToResMapFilePath = os.path.splitext(M2resFilePath)[0] + \
                    ".png"
                self.M2.showMirResMap(
                    resFile=M2resFilePath,
                    writeToResMapFilePath=writeToResMapFilePath)

        return pertCmdFilePath

    def _getPhoSimCamSurf(self, camSurfName):
        """Get the camera surface name used in PhoSim.

        Parameters
        ----------
        camSurfName : str
            Camera surface name (e.g. L1S1zer).

        Returns
        -------
        SurfaceType
            Camera surface enum.

        Raises
        ------
        ValueError
            Can not get the camera surface name.
        """

        # Get the camera surface name
        m = re.match(r"(\AL\d)S(\d)zer", camSurfName)
        if (m is None):
            raise ValueError("Cannot get the camera surface name: %s."
                             % camSurfName)

        # Get the surface name used in PhoSim
        camPhoFaceDict = {"1": "F", "2": "B"}
        surfName = m.groups()[0] + camPhoFaceDict[m.groups()[1]]

        return mapSurfNameToEnum(surfName)


if __name__ == "__main__":
    pass
