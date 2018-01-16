import os, re, shutil, unittest
import numpy as np

from wepPhoSim.CamSim import CamSim
from wepPhoSim.M2Sim import M2Sim
from wepPhoSim.M1M3Sim import M1M3Sim
from wepPhoSim.PhosimCommu import PhosimCommu

from wepPhoSim.OpdMetrology import OpdMetrology
from wepPhoSim.SkySim import SkySim

class TeleFacade(object):
    
    def __init__(self, cam=None, M1M3=None, M2=None, phoSimCommu=None, configFilePath=None):
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

    def getPhoSimArgs(self, instFilePath, cmdFilePath=None, numPro=1, numThread=1, outputDir=None, 
                        sensorName=None, e2ADC=1, logFilePath=None):
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

        argString = self.phoSimCommu.getPhoSimArgs(instFilePath, extraCommand=cmdFilePath, 
                                numProc=numPro, numThread=numThread, outputDir=outputDir, 
                                instrument=self.instName, sensorName=sensorName, e2ADC=e2ADC, 
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

    def setSubSysConfigFile(self, camDataDir=None, M1M3dataDir=None, M2dataDir=None, phosimDir=None):
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

        if (self.M2 is not None) and (M2dataDir is not None):
            self.M2.setMirrorDataDir(M2dataDir)

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

    def writeCmdFile(self, cmdFileDir, cmdSettingFile=None, pertFilePath=None, cmdFileName="taskPert.cmd"):
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
            self.phoSimCommu.writeToFile(cmdFilePath, sourceFile=cmdSettingFile)

        # Add the subsystem perturbation
        if (pertFilePath is not None):
            self.phoSimCommu.writeToFile(cmdFilePath, sourceFile=pertFilePath)

        return cmdFilePath

    def writeStarInstFile(self, instFileDir, skySim, obsId, aFilter, boresight=(0,0), camRot=0, 
                            mjd=59552.3, sedName="sed_500.txt", sciSensorOn=False, wfSensorOn=False, 
                            guidSensorOn=False, instSettingFile=None, instFileName="star.inst"):
        """
        
        Write the star instance file.
        
        Arguments:
            instFileDir {[str]} -- Directory to instance file.
            skySim {[SkySim]} -- SkySim object.
            obsId {[int]} -- Observation ID.
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
        
        Keyword Arguments:
            boresight {tuple} -- Telescope boresight in (ra, decl). (default: {(0,0)})
            camRot {float} -- Camera rotation angle. (default: {0})
            mjd {float} -- MJD of observation. (default: {59552.3})
            sedName {str} -- The name of the SED file with a file path that is relative to the 
                             data directory in PhoSim. (default: {"sed_500.txt"})
            sciSensorOn {bool} -- Science sensors are on. (default: {False})
            wfSensorOn {bool} -- Wavefront sensors are on. (default: {False})
            guidSensorOn {bool} -- Guider sensors are on. (default: {False})
            instSettingFile {[str]} -- Instance setting file. (default: {"None"}) 
            instFileName {str} -- Star instance file name. (default: {"star.inst"})
        
        Returns:
            [str] -- Instance file path.
        """

        # Instance file path
        instFilePath = os.path.join(instFileDir, instFileName)

        # Get the filter ID in PhoSim
        aFilterId = self.phoSimCommu.getFilterId(aFilter)

        # Write the default instance setting
        ra = boresight[0]
        dec = boresight[1]
        self.phoSimCommu.getStarInstance(obsId, aFilterId, ra=ra, dec=dec, rot=-camRot, 
                                            mjd=mjd, filePath=instFilePath)
        if (instFilePath is not None):
            self.phoSimCommu.writeToFile(instFilePath, sourceFile=instSettingFile)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Set the camera configuration
        content += self.phoSimCommu.doCameraConfig(sciSensorOn=sciSensorOn, wfSensorOn=wfSensorOn, 
                                                    guidSensorOn=guidSensorOn)

        # Write the star source
        for ii in range(len(skySim.starId)):
            content += self.phoSimCommu.generateStar(skySim.starId[ii], skySim.ra[ii], 
                                                skySim.decl[ii], skySim.mag[ii], sedName)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        return instFilePath

    def writeOpdInstFile(self, instFileDir, opdMetr, obsId, aFilter, wavelengthInNm, 
                            instSettingFile=None, instFileName="opd.inst"):
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
        aFilterId = self.phoSimCommu.getFilterId(aFilter)

        # Write the default instance setting
        self.phoSimCommu.getOpdInstance(obsId, aFilterId, filePath=instFilePath)
        if (instSettingFile is not None):
            self.phoSimCommu.writeToFile(instFilePath, sourceFile=instSettingFile)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Write the OPD source
        for ii in range(len(opdMetr.fieldX)):
            content += self.phoSimCommu.generateOpd(ii, opdMetr.fieldX[ii], opdMetr.fieldY[ii], wavelengthInNm)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        # Write the OPD SED file if necessary
        if (self.phoSimCommu.phosimDir is not None):
            self.phoSimCommu.writeSedFile(wavelengthInNm)
        else:
            print("Do not inspect the SED file for no setting of PhoSim directory.")

        return instFilePath

    def writePertBaseOnConfigFile(self, pertCmdFileDir, zAngleInDeg=0, rotAngInDeg=0, seedNum=None, 
                                    saveResMapFig=False, pertCmdFileName="pert.cmd"):
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
                randSurfInM = self.M1M3.genMirSurfRandErr(zAngleInRad, seedNum=seedNum)

            # Do the temperature correction
            M1M3TBulk = self.getConfigValue("M1M3TBulk")
            M1M3TxGrad = self.getConfigValue("M1M3TxGrad")
            M1M3TyGrad = self.getConfigValue("M1M3TyGrad")
            M1M3TzGrad = self.getConfigValue("M1M3TzGrad")
            M1M3TrGrad = self.getConfigValue("M1M3TrGrad")
            tempCorrInUm = self.M1M3.getTempCorr(M1M3TBulk, M1M3TxGrad, M1M3TyGrad, M1M3TzGrad, M1M3TrGrad)

            # Set the mirror surface in mm
            if (randSurfInM is not None):
                mirrorSurfInUm = (printthzInM + randSurfInM)*1e6 + tempCorrInUm
            else:
                mirrorSurfInUm = printthzInM*1e6 + tempCorrInUm
            self.M1M3.setSurfAlongZ(mirrorSurfInUm)

            resFile = [M1resFilePath, M3resFilePath]
            self.M1M3.writeMirZkAndGridResInZemax(resFile=resFile, numTerms=numTerms, writeZcInMnToFilePath=M1M3zcFilePath)

            # Get the Zk in mm
            zkInMm = np.loadtxt(M1M3zcFilePath)

            # Do the surface perturbation
            surfList = ["M1", "M3"]
            surfIdList = []
            for ii in range(2):
                surf = surfList[ii]
                surfId = self.phoSimCommu.getSurfaceId(surf)
                content += self.phoSimCommu.doSurfPert(surfId, zkInMm)

                # Collect the surface ID
                surfIdList.append(surfId)

                # Do the surface residue map perturbation
                content += self.phoSimCommu.doSurfMapPert(surfId, resFile[ii], 1)

            # Do the surface linkage
            content += self.phoSimCommu.doSurfLink(surfIdList[1], surfIdList[0])

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
            self.M2.writeMirZkAndGridResInZemax(resFile=M2resFilePath, surfaceGridN=surfaceGridN, 
                                            numTerms=numTerms, writeZcInMnToFilePath=M2zcFilePath)

            # Get the Zk in mm
            zkInMm = np.loadtxt(M2zcFilePath)

            # Do the surface perturbation
            surfId = self.phoSimCommu.getSurfaceId("M2")
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
            distTypeList = ["L1S1zer", "L1S2zer", "L2S1zer", "L2S2zer", "L3S1zer", "L3S2zer"]

            # Write the perturbation file
            for distType in distTypeList:
                # Get the surface ID
                surfName = self.__getPhoSimCamSurfName(distType)
                surfId = self.phoSimCommu.getSurfaceId(surfName)

                # Do the perturbation
                zkInMm = self.cam.getCamDistortionInMm(zAngleInRad, distType)
                content += self.phoSimCommu.doSurfPert(surfId, zkInMm)

        # Write the perturbation command to file
        self.phoSimCommu.writeToFile(pertCmdFilePath, content=content, mode="w")

        # Save the mirror residue map if necessary
        if (saveResMapFig):
            if (self.M1M3 is not None):
                resFile = [M1resFilePath, M3resFilePath]
                writeToResMapFilePath1 = os.path.splitext(M1resFilePath)[0] + ".png"
                writeToResMapFilePath3 = os.path.splitext(M3resFilePath)[0] + ".png"
                writeToResMapFilePath = [writeToResMapFilePath1, writeToResMapFilePath3]
                self.M1M3.showMirResMap(numTerms=numTerms, resFile=resFile, 
                                        writeToResMapFilePath=writeToResMapFilePath)

            if (self.M2 is not None):
                writeToResMapFilePath = os.path.splitext(M2resFilePath)[0] + ".png"
                self.M2.showMirResMap(numTerms=numTerms, resFile=M2resFilePath, 
                                      writeToResMapFilePath=writeToResMapFilePath)

        return pertCmdFilePath
                
    def __getPhoSimCamSurfName(self, camSurfName):
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

        return surfName

class TeleFacadeTest(unittest.TestCase):
    
    """
    Test functions in TeleFacade.
    """

    def setUp(self):

        # Set the configuration file path
        self.configFilePath = os.path.join("..", "data", "telescopeConfig", "GT.inst")

        # Set the subsystem data directory
        self.camDataDir = os.path.join("..", "data", "camera")
        self.M1M3dataDir = os.path.join("..", "data", "M1M3")
        self.M2dataDir = os.path.join("..", "data", "M2")
        self.phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"

        # Set the output dir
        self.outputDir = os.path.join("..", "output", "temp")
        os.makedirs(self.outputDir)

        # Set the command setting file
        self.starCmdSettingFile = os.path.join("..", "data", "cmdFile", "starDefault.cmd")

        # Set the instance setting file
        self.starInstSettingFile = os.path.join("..", "data", "instFile", "starDefault.inst")
        self.opdInstSettingFile = os.path.join("..", "data", "instFile", "opdDefault.inst")

    def testFunc(self):
        
        # Instantiate the needed objects
        cam = CamSim()
        M1M3 = M1M3Sim()
        M2 = M2Sim()
        phoSimCommu = PhosimCommu()
        metr = OpdMetrology()
        skySim = SkySim()

        # Instantiate the telescope facade class
        tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

        # Set the configuration file path
        tele.setConfigFile(self.configFilePath)
        self.assertEqual(tele.configFile, self.configFilePath)

        instFilePath = "temp.inst"
        argString = tele.getPhoSimArgs(instFilePath)
        ansArgString = "%s -i lsst -e 1" % os.path.abspath(instFilePath)
        self.assertEqual(argString, ansArgString)

        instruFile = "comcam10"
        tele.setInstName(instruFile)
        self.assertEqual(tele.instName, "comcam")
        self.assertEqual(tele.defocalDisInMm, 1.0)

        dofInUm = np.random.rand(50)
        tele.setDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(tele.dofInUm-dofInUm)), 0)

        tele.setDofInUm(np.zeros(50))
        tele.accDofInUm(dofInUm)
        self.assertEqual(np.sum(np.abs(tele.dofInUm-dofInUm)), 0)

        tele.setSubSysConfigFile(camDataDir=self.camDataDir, M1M3dataDir=self.M1M3dataDir, 
                                 M2dataDir=self.M2dataDir)
        self.assertEqual(tele.cam.camDataDir, self.camDataDir)
        self.assertEqual(tele.phoSimCommu.phosimDir, None)

        varName = "M1M3TxGrad"
        value = tele.getConfigValue(varName)
        self.assertEqual(value, -0.0894)

        dofFilePath = tele.writeAccDofFile(self.outputDir)
        self.assertLess(np.sum(np.abs(np.loadtxt(dofFilePath)-dofInUm)), 1e-7)
        os.remove(dofFilePath)

        zAngleInDeg = 27.0912
        rotAngInDeg = -1.2323/np.pi*180.0
        iSim = 6
        # pertCmdFilePath = "/Users/Wolf/Documents/stash/ts_tcs_wep_phosim/output/pert.cmd"
        pertCmdFilePath = tele.writePertBaseOnConfigFile(self.outputDir, zAngleInDeg=zAngleInDeg, 
                                           rotAngInDeg=rotAngInDeg, seedNum=iSim, saveResMapFig=True)
        cmdFile = open(pertCmdFilePath, "r")
        lines = cmdFile.readlines()
        cmdFile.close()
        self.assertEqual(len(lines), 256)

        cmdFilePath = tele.writeCmdFile(self.outputDir, cmdSettingFile=self.starCmdSettingFile, 
                                        pertFilePath=pertCmdFilePath, cmdFileName="star.cmd")

        starCmdFile = open(cmdFilePath, "r")
        lines = starCmdFile.readlines()
        starCmdFile.close()
        self.assertEqual(len(lines), 267)

        metr.addFieldXYbyDeg(0, 0)
        obsId = 9006000
        aFilter = "g"
        wavelengthInNm = 500
        instFilePath = tele.writeOpdInstFile(self.outputDir, metr, obsId, aFilter, wavelengthInNm, 
                                                instSettingFile=self.opdInstSettingFile)
        opdInstFile = open(instFilePath, "r")
        lines = opdInstFile.readlines()
        opdInstFile.close()
        self.assertEqual(len(lines), 55)

        skySim.addStarByRaDecInDeg(0, 1.0, 1.0, 17.0)
        boresight = (0.2, 0.3)
        instFilePath = tele.writeStarInstFile(self.outputDir, skySim, obsId, aFilter, boresight, 
                                                wfSensorOn=True, instSettingFile=self.starInstSettingFile)

        starInstFile = open(instFilePath, "r")
        lines = starInstFile.readlines()
        starInstFile.close()
        self.assertEqual(len(lines), 62)

        shutil.rmtree(self.outputDir)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
