import os, unittest, re
import numpy as np

from wepPhoSim.CamSim import CamSim
from wepPhoSim.M2Sim import M2Sim
from wepPhoSim.M1M3Sim import M1M3Sim
from wepPhoSim.PhosimCommu import PhosimCommu

class TeleFacade(object):
    
    def __init__(self, cam=None, M1M3=None, M2=None, phoSimCommu=None, configFilePath=None):
        
        self.cam = cam
        self.M1M3 = M1M3
        self.M2 = M2
        self.phoSimCommu = phoSimCommu

        self.dofInUm = np.zeros(50)

        self.configFile = configFilePath

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

    def setSubSysConfigFile(self, camDataDir=None, M1M3dataDir=None, M2dataDir=None):
        """
        
        Set the subsystem data directory.
        
        Keyword Arguments:
            camDataDir {[str]} -- Camera data directory. (default: {None})
            M1M3dataDir {[str]} -- M1M3 data directory. (default: {None})
            M2dataDir {[str]} -- M2 data directory. (default: {None})
        """
        
        if (self.cam is not None) and (camDataDir is not None):
            self.cam.setCamDataDir(camDataDir)

        if (self.M1M3 is not None) and (M1M3dataDir is not None):
            self.M1M3.setMirrorDataDir(M1M3dataDir)

        if (self.M2 is not None) and (M2dataDir is not None):
            self.M2.setMirrorDataDir(M2dataDir)

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

    def writeOpdInstFile(self, instFileDir):
        pass

    def writeStarInstFile(self, instFileDir):
        pass

    def writePertBaseOnConfigFile(self, pertCmdFileDir, seedNum=None, saveResMapFig=False):
        """
        
        Write the perturbation command file based on the telescope configuration file.
        
        Arguments:
            pertCmdFileDir {[str]} -- Directory to the pertubation command file. 
        
        Keyword Arguments:
            seedNum {[int]} -- Random seed number. (default: {None})
            saveResMapFig {[bool]} -- Save the mirror surface residue map or not. (default: {False}) 
        
        Returns:
            [str] -- Perturbation commend content.
        """

        # Perturbation command file name
        pertCmdFileName = "pert.cmd"

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
            printthzInM = M1M3.getPrintthz(zAngleInRad)

            # Add the surface error if necessary
            randSurfInM = None
            if (seedNum is not None):
                randSurfInM = M1M3.genMirSurfRandErr(zAngleInRad, seedNum=seedNum)

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

        return content
                
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

if __name__ == "__main__":
    
    cam = CamSim()
    M1M3 = M1M3Sim()
    M2 = M2Sim()
    phoSimCommu = PhosimCommu()

    tele = TeleFacade(cam=cam, M1M3=M1M3, M2=M2, phoSimCommu=phoSimCommu)

    # Set the telescope degree of freedom

    # Subsystem data direction
    camDataDir = "../data/camera"
    M1M3dataDir = "../data/M1M3"
    M2dataDir = "../data/M2"

    # Set the configuration file path
    configFilePath = "../data/telescopeConfig/GT.inst"
    tele.setConfigFile(configFilePath)

    # Get the value
    varName = "zenithAngle"
    value = tele.getConfigValue(varName, index=0)

    # Set the subsystem directory
    tele.setSubSysConfigFile(camDataDir=camDataDir, M1M3dataDir=M1M3dataDir, M2dataDir=M2dataDir)

    # Set the path of perturbation file
    iSim = 6
    pertFileDir = "../output"

    # Write the perturbation command file
    content = tele.writePertBaseOnConfigFile(pertFileDir, seedNum=iSim, saveResMapFig=True)

    # Write the opd instance command file


