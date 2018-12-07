import os, re
import numpy as np

from lsst.ts.phosim.Utility import SurfaceType, CamDistType


class TeleFacade(object):
    
    def __init__(self, cam=None, M1M3=None, M2=None, phoSimCommu=None,
                 configFilePath=None):
        """
        
        Initiate the TeleFacade object. This class uses the facade pattern that the high level
        class telescope helps to write the perturbations of camera, M1M3, and M2 into the PhoSim 
        by the interface to PhoSim. 
        
        Keyword Arguments:
            cam {[CamSim]} -- CamSim object. (default: {None})
            M1M3 {[M1M3Sim]} -- M1M3Sim object. (default: {None})
            M2 {[M2Sim]} -- M2Sim object. (default: {None})
            phoSimCommu {[PhosimCommu]} -- PhosimCommu object. (default: {None})
            configFilePath {[str]} -- Telescope configuration file path. (default: {None})
        """
        
        self.cam = cam
        self.M1M3 = M1M3
        self.M2 = M2
        self.phoSimCommu = phoSimCommu

        self.instName = "lsst"
        self.defocalDisInMm = 1.5
        self.dofInUm = np.zeros(50)

        self.configFile = configFilePath

    def runPhoSim(self, argString):
        """
        
        Run the PhoSim program.
        
        Arguments:
            argString {[str]} -- Arguments for PhoSim.
        """
        
        self.phoSimCommu.runPhoSim(argstring=argString)

    def getPhoSimArgs(self, instFilePath, cmdFilePath=None, numPro=1,
                      numThread=1, outputDir=None, sensorName=None,
                      e2ADC=1, logFilePath=None):
        """
        
        Get the arguments needed to run the PhoSim.
        
        Arguments:
            instFilePath {[str]} -- Instance catalog file path.
        
        Keyword Arguments:
            cmdFilePath {[str]} -- Command file to modify the default physics. (default: {None})
            numPro {int} -- Number of processors. (default: {1})
            numThread {int} -- Number of threads. (default: {1})
            outputDir {[str]} -- Output image directory. (default: {None})
            sensorName {str} -- Sensor chip specification (e.g., all, R22_S11, "R22_S11|R22_S12") 
                                (default: {None})
            e2ADC {int} -- Whether to generate amplifier images (1 = true, 0 = false). (default: {1})
            logFilePath {[str]} -- Log file path for PhoSim calculation log. (default: {None})
        
        Returns:
            [str] -- Arguments to run the PhoSim.
        """

        argString = self.phoSimCommu.getPhoSimArgs(
                            instFilePath, extraCommandFile=cmdFilePath, 
                            numProc=numPro, numThread=numThread,
                            outputDir=outputDir, instrument=self.instName,
                            sensorName=sensorName, e2ADC=e2ADC, 
                            logFilePath=logFilePath)

        return argString

    def setInstName(self, instruFile):
        """
        
        Set the instrument name from the instrument file.
        
        Arguments:
            instruFile {[str]} -- Instrument folder name.
        
        Returns:
            [str] -- Instrument name.
            [float] -- Defocal offset in mm.
        
        Raises:
            RuntimeError -- No instrument found.
        """

        # Get the instrument name
        m = re.match(r"([a-z]+)(?:(\d+))?$", instruFile)
        if m is None:
             raise RuntimeError("Cannot get the instrument name: %s." % instruFile)
        instName = m.groups()[0]

        # Decide the defocal distance offset in mm
        defocalOffset = m.groups()[1]
        if (defocalOffset is not None):
            defocalOffset = float(defocalOffset)/10
        else:
            # Default defocal distance is 1.5 mm
            defocalOffset = 1.5

        self.instName = instName
        self.defocalDisInMm = defocalOffset

    def setDofInUm(self, dofInUm):
        """
        
        Set the accumulated degree of freedom (DOF) in um.
        idx 0-4: M2 dz, dx, dy, rx, ry
        idx 5-9: Cam dz, dx, dy, rx, ry
        idx 10-29: M1M3 20 bending modes
        idx 30-49: M2 20 bending modes
        
        Arguments:
            dofInUm {[list/ ndarray]} -- DOF in um.
        """

        self.dofInUm = dofInUm

    def accDofInUm(self, dofInUm):
        """
        
        Accumulate the aggregated degree of freedom (DOF) in um.
        idx 0-4: M2 dz, dx, dy, rx, ry
        idx 5-9: Cam dz, dx, dy, rx, ry
        idx 10-29: M1M3 20 bending modes
        idx 30-49: M2 20 bending modes
        
        Arguments:
            dofInUm {[list/ ndarray]} -- DOF in um.
        """

        self.dofInUm += dofInUm

    def setSubSysConfigFile(self, camDataDir=None, M1M3dataDir=None,
                            M2dataDir=None, phosimDir=None):
        """
        
        Set the subsystem data directory.
        
        Keyword Arguments:
            camDataDir {[str]} -- Camera data directory. (default: {None})
            M1M3dataDir {[str]} -- M1M3 data directory. (default: {None})
            M2dataDir {[str]} -- M2 data directory. (default: {None})
            phosimDir {[str]} -- PhoSim directory. (default: {None})
        """
        
        if (self.cam is not None) and (camDataDir is not None):
            self.cam.setCamDataDir(camDataDir)

        if (self.M1M3 is not None) and (M1M3dataDir is not None):
            self.M1M3.setMirrorDataDir(M1M3dataDir)
            self.M1M3.config()

        if (self.M2 is not None) and (M2dataDir is not None):
            self.M2.setMirrorDataDir(M2dataDir)
            self.M2.config()

        if (self.phoSimCommu is not None) and (phosimDir is not None):
            self.phoSimCommu.setPhoSimDir(phosimDir)

    def setConfigFile(self, configFilePath):
        """
        
        Set the telescope configuration file.
        
        Arguments:
            configFilePath {[str]} -- Configuration file path.
        """

        self.configFile = configFilePath

    def getConfigValue(self, varName, index=1):
        """
        
        Get the value of certain variable defined in the configuration file.
        
        Arguments:
            varName {[str]} -- Name of variable.
        
        Keyword Arguments:
            index {int} -- Index of value. (default: {1})
        
        Returns:
            [float/ int/ str] -- Variable value.
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

    def writeAccDofFile(self, outputFileDir, dofFileName="pert.mat"):
        """
        
        Write the accumulated degree of freedom (DOF) in um to file.
        
        Arguments:
            outputFileDir {[str]} -- Output file directory.

        Keyword Arguments:
            dofFileName {[str]} -- DOF file name. (default: {"pert.mat"})

        Returns:
            [str] -- DOF file path.
        """
        
        dofFilePath = os.path.join(outputFileDir, dofFileName)
        np.savetxt(dofFilePath, self.dofInUm)

        return dofFilePath

    def writeCmdFile(self, cmdFileDir, cmdSettingFile=None, pertFilePath=None,
                     cmdFileName="taskPert.cmd"):
        """
        
        Write the physical command file.
        
        Arguments:
            cmdFileDir {[str]} -- Directory to the OPD command file.
        
        Keyword Arguments:
            cmdSettingFile {[str]} -- Physical command setting file path. (default: {None})
            pertFilePath {[str]} -- Subsystem perturbation command file path. (default: {None})
            cmdFileName {str} -- Command file name. (default: {"taskPert.cmd"})
        
        Returns:
            [str] -- Command file path.
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

    def writeStarInstFile(self, instFileDir, skySim, obsId, filterType,
                          boresight=(0,0), rot=0, mjd=59552.3, simSeed=1000, 
                          sedName="sed_500.txt", sciSensorOn=False, 
                          wfSensorOn=False, guidSensorOn=False,
                          instSettingFile=None, instFileName="star.inst"):
        """
        
        Write the star instance file.
        
        Arguments:
            instFileDir {[str]} -- Directory to instance file.
            skySim {[SkySim]} -- SkySim object.
            obsId {[int]} -- Observation ID.
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
        
        Keyword Arguments:
            boresight {[tuple]} -- Telescope boresight in (ra, decl). (default: {(0,0)})
            rot {[float]} -- Angle of sky relative to camera coordinates (from North over East) in 
                           decimal degrees. (default: {0})
            mjd {[float]} -- MJD of observation. (default: {59552.3})
            simSeed {[int]} -- Random number seed. (default: {1000})
            sedName {[str]} -- The name of the SED file with a file path that is relative to the 
                             data directory in PhoSim. (default: {"sed_500.txt"})
            sciSensorOn {[bool]} -- Science sensors are on. (default: {False})
            wfSensorOn {[bool]} -- Wavefront sensors are on. (default: {False})
            guidSensorOn {[bool]} -- Guider sensors are on. (default: {False})
            instSettingFile {[str]} -- Instance setting file. (default: {"None"}) 
            instFileName {[str]} -- Star instance file name. (default: {"star.inst"})
        
        Returns:
            [str] -- Instance file path.
        """

        # Instance file path
        instFilePath = os.path.join(instFileDir, instFileName)

        # Get the filter ID in PhoSim
        aFilterId = self.phoSimCommu.getFilterId(filterType)

        # Write the default instance setting
        ra = boresight[0]
        dec = boresight[1]
        self.phoSimCommu.getStarInstance(
                            obsId, aFilterId, ra=ra, dec=dec, rot=rot, 
                            mjd=mjd, simSeed=simSeed, filePath=instFilePath)
        if (instFilePath is not None):
            self.phoSimCommu.writeToFile(instFilePath,
                                         sourceFile=instSettingFile)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Set the camera configuration
        content += self.phoSimCommu.doCameraConfig(
                                        sciSensorOn=sciSensorOn,
                                        wfSensorOn=wfSensorOn,
                                        guidSensorOn=guidSensorOn)

        # Write the star source
        for ii in range(len(skySim.starId)):
            content += self.phoSimCommu.generateStar(
                                            skySim.starId[ii], skySim.ra[ii], 
                                            skySim.decl[ii], skySim.mag[ii],
                                            sedName)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        return instFilePath

    def writeOpdInstFile(self, instFileDir, opdMetr, obsId, filterType,
                         wavelengthInNm, instSettingFile=None,
                         instFileName="opd.inst"):
        """
        
        Write the optical path difference (OPD) instance file.
        
        Arguments:
            instFileDir {[str]} -- Directory to instance file.
            opdMetr {[OpdMetrology]} -- OpdMetrology object.
            obsId {[int]} -- Observation ID.
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
            wavelengthInNm {[float]} -- OPD source wavelength in nm.

        Keyword Arguments:
            instSettingFile {[str]} -- Instance setting file. (default: {"None"}) 
            instFileName {[str]} -- OPD instance file name. (default: {"opd.inst"})
        
        Returns:
            [str] -- Instance file path.
        """
        
        # Instance file path
        instFilePath = os.path.join(instFileDir, instFileName)

        # Get the filter ID in PhoSim
        aFilterId = self.phoSimCommu.getFilterId(filterType)

        # Write the default instance setting
        self.phoSimCommu.getOpdInstance(obsId, aFilterId,
                                        filePath=instFilePath)
        if (instSettingFile is not None):
            self.phoSimCommu.writeToFile(instFilePath,
                                         sourceFile=instSettingFile)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Write the OPD source
        for ii in range(len(opdMetr.fieldX)):
            content += self.phoSimCommu.generateOpd(
                                            ii, opdMetr.fieldX[ii],
                                            opdMetr.fieldY[ii],
                                            wavelengthInNm)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        # Write the OPD SED file if necessary
        if os.path.isdir(self.phoSimCommu.phosimDir):
            self.phoSimCommu.writeSedFile(wavelengthInNm)
        else:
            print("Do not inspect the SED file for no setting of PhoSim directory.")

        return instFilePath

    def writePertBaseOnConfigFile(self, pertCmdFileDir, zAngleInDeg=0,
                                  rotAngInDeg=0, seedNum=None, 
                                  saveResMapFig=False,
                                  pertCmdFileName="pert.cmd"):
        """
        
        Write the perturbation command file based on the telescope configuration file.
        
        Arguments:
            pertCmdFileDir {[str]} -- Directory to the pertubation command file. 
        
        Keyword Arguments:
            zAngleInDeg {[float]} -- Zenith angle in degree. (default: {0})
            rotAngInDeg {[float]} -- Camera rotation angle in degree between -90 and 90 degrees. 
                                     (default: {0})
            seedNum {[int]} -- Random seed number. (default: {None})
            saveResMapFig {[bool]} -- Save the mirror surface residue map or not. (default: {False})
            pertCmdFileName {[str]} -- Perturbation command file name. (default: {pert.cmd})
        
        Returns:
            [str] -- Perturbation command file path.
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
        zAngleInRad = zAngleInDeg/180.0*np.pi

        # Get the number of zernike terms used in Zemax
        numTerms = self.getConfigValue("znPert")

        # Get the numeber of grid used in Zemax
        surfaceGridN = self.getConfigValue("surfaceGridN")

        # Write the camera perturbation command file
        if (self.M1M3 is not None):
            
            # Do the gravity correction
            printthzInM = self.M1M3.getPrintthz(zAngleInRad)

            # Add the surface error if necessary
            randSurfInM = None
            if (seedNum is not None):
                randSurfInM = self.M1M3.genMirSurfRandErr(zAngleInRad,
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
                mirrorSurfInUm = (printthzInM + randSurfInM)*1e6 + \
                                 tempCorrInUm
            else:
                mirrorSurfInUm = printthzInM*1e6 + tempCorrInUm
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
            for ii in range(2):
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
            self.cam.setRotAngInDeg(rotAngInDeg)

            # Set the temperature information
            tempInDegC = self.getConfigValue("camTB")
            self.cam.setBodyTempInDegC(tempInDegC)

            # Get the camera distortion
            # distTypeList = ["L1S1zer", "L1S2zer", "L2S1zer", "L2S2zer", "L3S1zer", "L3S2zer"]

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
                writeToResMapFilePath1 = os.path.splitext(M1resFilePath)[0] + ".png"
                writeToResMapFilePath3 = os.path.splitext(M3resFilePath)[0] + ".png"
                writeToResMapFilePath = [writeToResMapFilePath1, writeToResMapFilePath3]
                self.M1M3.showMirResMap(
                                resFile=resFile, 
                                writeToResMapFilePath=writeToResMapFilePath)

            if (self.M2 is not None):
                writeToResMapFilePath = os.path.splitext(M2resFilePath)[0] + ".png"
                self.M2.showMirResMap(
                                resFile=M2resFilePath, 
                                writeToResMapFilePath=writeToResMapFilePath)

        return pertCmdFilePath
                
    def _getPhoSimCamSurf(self, camSurfName):
        """
        
        Get the camera surface name used in PhoSim.
        
        Arguments:
            camSurfName {[str]} -- Camera surface name.
        
        Returns:
            [str] -- Camera surface name in PhoSim.
        
        Raises:
            RuntimeError -- Can not get the camera surface name.
        """

        # Get the camera surface name
        m = re.match(r"(\AL\d)S(\d)zer", camSurfName)
        if m is None:
            raise RuntimeError("Cannot get the camera surface name: %s." % camSurfName)
        
        # Get the surface name used in PhoSim
        camPhoFaceDict = {"1": "F", "2": "B"}
        surfName = m.groups()[0] + camPhoFaceDict[m.groups()[1]]

        # Modify this to get the surface type
        if (surfName == "L1F"):
            surfaceType = SurfaceType.L1F
        if (surfName == "L1B"):
            surfaceType = SurfaceType.L1B
        if (surfName == "L2F"):
            surfaceType = SurfaceType.L2F
        if (surfName == "L2B"):
            surfaceType = SurfaceType.L2B
        if (surfName == "L3F"):
            surfaceType = SurfaceType.L3F
        if (surfName == "L3B"):
            surfaceType = SurfaceType.L3B

        return surfaceType


if __name__ == "__main__":
    pass
