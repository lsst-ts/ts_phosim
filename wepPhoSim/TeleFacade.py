import os, re, unittest
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
            e2ADC {int} -- Whether to generate amplifier images (1 = true, 0 = false). (default: {1})
            logFilePath {[str]} -- Log file path for PhoSim calculation log. (default: {None})
        
        Returns:
            [str] -- Arguments to run the PhoSim.
        """

        argString = self.phoSimCommu.getPhoSimArgs(instFilePath, extraCommand=cmdFilePath, 
                                numProc=numPro, numThread=numThread, outputDir=outputDir, 
                                instrument=self.instName, e2ADC=e2ADC, logFilePath=logFilePath)

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
        
        Arguments:
            dofInUm {[list/ ndarray]} -- DOF in um.
        """

        self.dofInUm = dofInUm

    def accDofInUm(self, dofInUm):
        """
        
        Accumulate the aggregated degree of freedom (DOF) in um.
        
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

    def writeDefaultStarCmdFile(self, cmdFileDir, pertFilePath=None, cmdFileName="starPert.cmd"):
        """
        
        Write the default physical command file of stars with the subsystem perturbation 
        if necessary. 
        
        Arguments:
            cmdFileDir {[str]} -- Directory to the star command file.
        
        Keyword Arguments:
            cmdFileName {[str]} -- Star command file name. (default: {"starPert.cmd"})
            pertFilePath {[str]} -- Subsystem perturbation command file path. (default: {None})
        
        Returns:
            [str] -- Command file path.
        """
        
        # Command file path
        cmdFilePath = os.path.join(cmdFileDir, cmdFileName)

        # Write the default physical setting
        self.phoSimCommu.getDefaultCmd(filePath=cmdFilePath)

        # Add the subsystem perturbation
        if (pertFilePath is not None):
            self.phoSimCommu.writeToFile(cmdFilePath, sourceFile=pertFilePath)

        return cmdFilePath

    def writeDefaultOpdCmdFile(self, cmdFileDir, pertFilePath=None, cmdFileName="opdPert.cmd"):
        """
        
        Write the default physical command file of optical path difference (OPD) with the 
        subsystem perturbation if necessary. 
        
        Arguments:
            cmdFileDir {[str]} -- Directory to the OPD command file.
        
        Keyword Arguments:
            cmdFileName {[str]} -- OPD command file name. (default: {"opdPert.cmd"})
            pertFilePath {[str]} -- Subsystem perturbation command file path. (default: {None})

        Returns:
            [str] -- Command file path.
        """

        # Command file path
        cmdFilePath = os.path.join(cmdFileDir, cmdFileName)

        # Write the default physical setting
        self.phoSimCommu.getDefaultOpdCmd(filePath=cmdFilePath)

        # Add the subsystem perturbation
        if (pertFilePath is not None):
            self.phoSimCommu.writeToFile(cmdFilePath, sourceFile=pertFilePath)

        return cmdFilePath

    def writeDefaultStarInstFile(self, instFileDir, skySim, obsId, aFilter, boresight=(0,0), 
                                    camRot=0, mjd=59552.3, sedName="sed_500.txt", sciSensorOn=False, 
                                    wfSensorOn=False, guidSensorOn=False, instFileName="star.inst"):
        """
        
        Write the default star instance file.
        
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
        self.phoSimCommu.getDefaultInstance(obsId, aFilterId, ra=ra, dec=dec, rot=-camRot, 
                                            mjd=mjd, filePath=instFilePath)

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

    def writeDefaultOpdInstFile(self, instFileDir, opdMetr, obsId, aFilter, wavelengthInNm, 
                                instFileName="opd.inst"):
        """
        
        Write the default optical path difference (OPD) instance file.
        
        Arguments:
            instFileDir {[str]} -- Directory to instance file.
            opdMetr {[OpdMetrology]} -- OpdMetrology object.
            obsId {[int]} -- Observation ID.
            aFilter {[str]} -- Active filter type ("u", "g", "r", "i", "z", "y").
            wavelengthInNm {[float]} -- OPD source wavelength in nm.

        Keyword Arguments:
            instFileName {[str]} -- OPD instance file name. (default: {"opd.inst"})
        
        Returns:
            [str] -- Instance file path.
        """
        
        # Instance file path
        instFilePath = os.path.join(instFileDir, instFileName)

        # Get the filter ID in PhoSim
        aFilterId = self.phoSimCommu.getFilterId(aFilter)

        # Write the default instance setting
        self.phoSimCommu.getDefaultOpdInstance(obsId, aFilterId, filePath=instFilePath)

        # Write the telescope accumulated degree of freedom (DOF)
        content = self.phoSimCommu.doDofPert(self.dofInUm)

        # Write the OPD source
        for ii in range(len(opdMetr.fieldX)):
            content += self.phoSimCommu.generateOpd(ii, opdMetr.fieldX[ii], opdMetr.fieldY[ii], wavelengthInNm)
        self.phoSimCommu.writeToFile(instFilePath, content=content)

        # Write the OPD SED file if necessary
        self.phoSimCommu.writeSedFile(wavelengthInNm)

        return instFilePath

    def writeAccDofFile(self, outputFileDir):
        """
        
        Write the accumulated degree of freedom (DOF) in um to file.
        
        Arguments:
            outputFileDir {[str]} -- Output file directory.
        """
        
        dofFileName = "pert.mat"
        dofFilePath = os.path.join(outputFileDir, dofFileName)
        np.savetxt(dofFilePath, self.dofInUm)

    def writePertBaseOnConfigFile(self, pertCmdFileDir, seedNum=None, saveResMapFig=False, pertCmdFileName="pert.cmd"):
        """
        
        Write the perturbation command file based on the telescope configuration file.
        
        Arguments:
            pertCmdFileDir {[str]} -- Directory to the pertubation command file. 
        
        Keyword Arguments:
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

        # Get the zenith angle in degree
        zAngleInDeg = self.getConfigValue("zenithAngle")
        zAngleInRad = zAngleInDeg/180.0*np.pi

        # Get the camera rotation angle
        rotAngInRad = self.getConfigValue("camRotation")

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
            self.cam.setRotAngInRad(rotAngInRad)

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
        pass

    def testFunc(self):
        pass



if __name__ == "__main__":

    # Do the unit test
    unittest.main()
    
    # cam = CamSim()
    # M1M3 = M1M3Sim()
    # M2 = M2Sim()
    # phoSimCommu = PhosimCommu()
    # metr = OpdMetrology()

    # tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # # Set the telescope degree of freedom

    # # Subsystem data direction
    # camDataDir = "../data/camera"
    # M1M3dataDir = "../data/M1M3"
    # M2dataDir = "../data/M2"
    # phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"

    # # Set the configuration file path
    # configFilePath = "../data/telescopeConfig/GT.inst"
    # tele.setConfigFile(configFilePath)

    # # Get the value
    # varName = "zenithAngle"
    # value = tele.getConfigValue(varName, index=1)

    # # Set the subsystem directory
    # tele.setSubSysConfigFile(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir, 
    #                          phosimDir=phosimDir)

    # # Set the path of perturbation file
    # iSim = 6
    # outputFileDir = "../output"

    # # Write the perturbation command file
    # # content = tele.writePertBaseOnConfigFile(outputFileDir, seedNum=iSim, saveResMapFig=True)[0]

    # # Set the default LSST GQ metrology
    # metr.setDefaultLsstGQ()

    # # Add the wavefront sensor
    # fieldWFSx, fieldWFSy = metr.getDefaultLsstWfsGQ()
    # metr.addFieldXYbyDeg(fieldWFSx, fieldWFSy)

    # # Write the opd physical command file
    # pertFilePath = os.path.join(outputFileDir, "pert.cmd")
    # cmdFilePath = tele.writeDefaultOpdCmdFile(outputFileDir, pertFilePath=pertFilePath)

    # # Write the opd instance file
    # obsId = 9006000
    # aFilter = "g"
    # wavelengthInNm = 500
    # instFilePath = tele.writeDefaultOpdInstFile(outputFileDir, metr, obsId, aFilter, wavelengthInNm)

    # # Write the accumulated DOF file
    # tele.writeAccDofFile(outputFileDir)

    # # Get the argument to run the phosim
    # outputImgFileDir = "../outputImg"
    # logFilePath = os.path.join(outputImgFileDir, "phosimOpd.log")
    # argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgFileDir, 
    #                                 e2ADC=0, logFilePath=logFilePath)
    
    # # Run the phosim
    # # tele.runPhoSim(argString)

    # # Calculate the PSSN
    # # zen = tele.getConfigValue("zenithAngle", index=1)
    # # opdFitsFile = os.path.join(outputImgFileDir, "opd_9006000_0.fits.gz")
    # # wavelengthInUm = 0.5
    # # pssn = metr.calcPSSN(opdFitsFile, wavelengthInUm, zen=zen, debugLevel=0)

    # # Write the default star command file
    # cmdFilePath = tele.writeDefaultStarCmdFile(outputFileDir, pertFilePath=pertFilePath)

    # # Set the sky information
    # skySim = SkySim()
    # skySim.addStarByRaDecInDeg(0, 1.196, 1.176, 17)

    # # Write the default star instance file
    # instFilePath = tele.writeDefaultStarInstFile(outputFileDir, skySim, obsId, aFilter, wfSensorOn=True, 
    #                                                 instFileName="star.inst")

    # # Get the argument to run the phosim
    # logFilePath = os.path.join(outputImgFileDir, "phosimStar.log")
    # argString = tele.getPhoSimArgs(instFilePath, cmdFilePath=cmdFilePath, numPro=2, outputDir=outputImgFileDir, 
    #                             e2ADC=0, logFilePath=logFilePath)

    # # Run the phosim
    # # tele.runPhoSim(argString)
