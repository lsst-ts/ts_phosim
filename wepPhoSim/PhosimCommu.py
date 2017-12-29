import os, subprocess
import numpy as np

class PhosimCommu(object):

    def __init__(self, phosimDir=None):
        """
        
        Initiate the object.
        
        Keyword Arguments:
            phosimDir {[str]} -- PhoSim directory. (default: {None})
        """
        
        self.phosimDir = phosimDir

    def setPhoSimDir(phosimDir):
        """
        
        Set the directory of PhoSim.
        
        Arguments:
            phosimDir {[str]} -- Directory of PhoSim.
        """

        self.phosimDir = phosimDir

    def getFilterId(self, aFilter):
        """
        
        Get the active filter ID in PhoSim.
        
        Arguments:
            aFilter {[str]} -- Active filter ("u", "g", "r", "i", "z", "y").
        
        Returns:
            [int] -- Active filter ID in PhoSim.
        """

        # Active filter ID in PhoSim
        phosimFilterID = {"u": 0, "g": 1, "r": 2, "i": 3, "z": 4, "y": 5}

        # Get the filter ID in PhoSim
        return self.__getMapId(phosimFilterID, aFilter)

    def getSurfaceId(self, aSurface):
        """
        
        Get the surface ID in PhoSim.
        
        Arguments:
            aSurface {[str]} -- Surface ("M1", "M2", "M3", "L1F", "L1B", "L2F", "L2B", 
                                "FilterF", "FilterB", "L3F", "L3B", "FP", "Chip")
        
        Returns:
            [int] -- Surface ID in PhoSim.
        """

        # LSST surface ID in PhoSim
        # "L1F" means lense 1 front. "L1B" means lense 1 back.
        # "chip" means individual chip.
        surfaceID = {"M1": 0, "M2": 1, "M3": 2, "L1F": 3, "L1B": 4, "L2F": 5, 
                     "L2B": 6, "FilterF": 7, "FilterB": 8, "L3F": 9, "L3B": 10,
                     "FP": 11, "Chip": 12}

        # Get the surface ID in PhoSim
        return self.__getMapId(surfaceID, aSurface)

    def __getMapId(self, dictMap, keyWord):
        """
        
        Get the ID in map dictionary defined in PhoSim.
        
        Arguments:
            dictMap {[dict]} -- Dictionary map.
            keyWord {[str]} -- Keyword.
        
        Returns:
            [int] -- ID.
        """

        # Get the key word ID in PhoSim
        keyId = None
        try:
            keyId = dictMap[keyWord]
        except KeyError:
            print("The keyWord '%s' is not supported." % keyWord)

        return keyId

    def doDofPert(self, dofInUm, phosimStrIdx=5, nDofInPhoSim=50):
        """
        
        Do the perturbation of degree of freedom (DOF).
        
        Arguments:
            dofInUm {[list/ ndarray]} -- DOF in um.
        
        Keyword Arguments:
            phosimStrIdx {int} -- Start index of DOF in PhoSim. (default: {5})
            nDofInPhoSim {int} -- Allowed number of DOF in PhoSim. (default: {50})
        
        Returns:
            [str] -- Perturbation command used in PhoSim.
        
        Raises:
            ValueError -- The number of DOF is incorrect.
        """
    
        # Check the length of degree of freedom.
        if (len(dofInUm) != int(nDofInPhoSim)):
            raise ValueError("The length of degree of freedom should be %s." % (int(nDofInPhoSim)))

        # List the phosim dof here
        # idx 5-9: M2 dz, dx, dy, rx, ry
        # idx 10-14: Cam dz, dx, dy, rx, ry
        # idx 15-24: M1M3 20 bending modes
        # idx 25-54: M2 20 bending modes

        # Write the perturbation of degree of freedom 
        content = ""
        for ii in range(len(dofInUm)):
            content += "move %d %7.4f \n" % (phosimStrIdx+ii, dofInUm[ii])

        return content

    def doSurfPert(self, surfId, zkInMm):
        """
        
        Do the perturbation of surface in the basis of Zernike polynomial.
        
        Arguments:
            surfId {[int]} -- PhoSim surface ID.
            zkInMm {[list/ ndarray]} -- Zernike terms in mm starts from z0.
        
        Returns:
            [str] -- Perturbation command used in PhoSim.
        """
        
        # Write the perturbation of surface
        content = ""
        for ii in range(len(zkInMm)):
            content += "izernike %d %d %s \n" % (surfId, ii, zkInMm[ii])

        return content

    def doSurfMapPert(self, surfId, surfFilePath, relScale):
        """
        
        Do the perturbation of surface map of individual optic.
        
        Arguments:
            surfId {[int]} -- PhoSim surface Id.
            surfFilePath {[str]} -- File path of surface.
            relScale {[int]} -- Relative scaling.
        
        Returns:
            [str] -- Perturbation command used in PhoSim.
        """

        # Write the perturbation of surface map
        content = "surfacemap %d %s %d \n" % (surfId, surfFilePath, relScale)

        return content

    def doCameraConfig(self, camConfigId):
        """
        
        Do the camera configuratioin.
        
        Arguments:
            camConfigId {[int]} -- Camera configuration. Bit mask that defines which sensor groups 
                                   are used. "1": science sensors on. "2": wavefront sensors on.
                                   "3": guiders on. "7": all sensors (default). 
        
        Returns:
            [str] -- Instance command used in PhoSim.
        """
        
        # Write the camera configuration
        content = "camconfig %d \n" % camConfigId

        return content

    def generateOpd(self, opdId, fieldXInDeg, fieldYInDeg, wavelengthInNm):
        """
        
        Generate the oprical path difference (OPD) source.
        
        Arguments:
            opdId {[int]} -- OPD ID.
            fieldXInDeg {[float]} -- Field X in degree.
            fieldYInDeg {[float]} -- Field Y in degree.
            wavelengthInNm {[float]} -- Wavelength of the OPD source in nm.
        
        Returns:
            [str] -- Perturbation command used in PhoSim.
        """

        # Write the OPD information
        content = "opd %2d\t%9.6f\t%9.6f %5.1f \n" % (opdId, fieldXInDeg, fieldYInDeg, wavelengthInNm)

        return content

    def generateStar(self, starId, ra, dec, magNorm, sedName, redshift=0, 
                     gamma1=0, gamma2=0, kappa=0, deltaRa=0, deltaDec=0, 
                     sourceType="star", spatialPars=0):
        """
        
        Generate the star source.
        
        Arguments:
            starId {[int]} -- Star ID.
            ra {[float]} -- The right ascension of the center of the object or image in decimal degrees.
            dec {[float]} -- The declination of the center of the object in decimal degrees.
            magNorm {[float]} -- The normalization of the flux of the object in AB magnitudes 
                                 at (500 nm)/(1+z) (which is roughly equivalent to V (AB) or g (AB)).
            sedName {[str]} -- The name of the SED file with a file path that is relative to the 
                               data directory in PhoSim.
        
        Keyword Arguments:
            redshift {number} -- The redshift (or blueshift) of the object. Note that the SED does not 
                                 need to be redshifted if using this. (default: {0})
            gamma1 {number} -- The value of the shear parameter gamma1 used in weak lensing. (default: {0})
            gamma2 {number} -- The value of the shear parameter gamma2 used in weak lensing. (default: {0})
            kappa {number} -- The value of the magnification parameter in weak lensing. (default: {0})
            deltaRa {number} -- The value of the declination offset in radians. This can be used either 
                                for weak lensing or objects that moved from another exposure if you do 
                                not want to change the source position in the first two columns. (default: {0})
            deltaDec {number} -- The value of the declination offset in radians. This can be used either for 
                                 weak lensing or objects that moved from another exposure if you do not want 
                                 to change the source position in the first two columns. (default: {0})
            sourceType {str} -- The name of the spatial model to be used as defined below. (default: {"star"})
            spatialPars {number} -- The associated parameters for each spatial model.There could be none or many. 
                                    While the parser is reading the model it looks for more parameters based on 
                                    the name of the model. (default: {0})
        
        Returns:
            [str] -- Perturbation command used in PhoSim.
        """

        # Write the star information
        content = "object %2d\t%9.6f\t%9.6f %9.6f ../sky/%s %.1f %.1f %.1f %.1f %.1f %.1f %s %.1f none none \n" % (
                    starId, ra, dec, magNorm, sedName, redshift, gamma1, gamma2, kappa, deltaRa, deltaDec, sourceType, 
                    spatialPars)

        return content

    def getDefaultCmd(self, filePath=None):
        """
        
        Get the default physical command.
        
        Keyword Arguments:
            filePath {[str]} -- File path to save the physical command. (default: {None})
        
        Returns:
            [str] -- Physical command used in PhoSim.
        """

        # Physics Commands
        content = ""
        content += "backgroundmode 0 \n"
        content += "raydensity 0.0 \n"
        content += "perturbationmode 1 \n"
        content += "trackingmode 0 \n"
        content += "cleartracking \n"
        content += "clearclouds \n"
        content += "lascatprob 0.0 \n"
        content += "contaminationmode 0 \n"
        content += "diffractionmode 1 \n"
        content += "straylight 0 \n"
        content += "detectormode 0 \n"

        if (filePath is not None):
            self.writeToFile(filePath, content=content, mode="w")

        return content

    def getDefaultOpdCmd(self, filePath=None):
        """
        
        Get the default optical path difference (OPD) physical command.
        
        Keyword Arguments:
            filePath {[str]} -- File path to save the physical command. (default: {None})
        
        Returns:
            [str] -- Physical command used in PhoSim.
        """
        
        # Physics Commands
        content = ""
        content += "backgroundmode 0 \n"
        content += "raydensity 0.0 \n"
        content += "perturbationmode 1 \n"

        if (filePath is not None):
            self.writeToFile(filePath, content=content, mode="w")

        return content

    def getDefaultInstance(self, obsId, aFilterId, ra=0, dec=0, rot=0, mjd=49552.3, filePath=None):
        """
        
        Get the default instance catalog.
        
        Arguments:
            obsId {[int]} -- Observation ID.
            aFilterId {[int]} -- PhoSim lsst filter ID.
        
        Keyword Arguments:
            ra {float} -- Unrefracted Right Ascension in decimal degrees. (default: {0})
            dec {float} -- Unrefracted Declination in decimal degrees. (default: {0})
            rot {float} -- Angle of sky relative to camera coordinates in ? (from North over East). 
                           (default: {0})
            mjd {float} -- MJD of observation. (default: {49552.3})
            filePath {[str]} -- File path to save the instance. (default: {None})
        
        Returns:
            [str] -- Instance command used in PhoSim.
        """

        # Observation Parameters
        content = ""
        content += "Opsim_obshistid %d \n" % obsId
        content += "Opsim_filter %d \n" % aFilterId
        content += "mjd %.10f \n" % mjd
        content += "SIM_VISTIME 33.0 \n"
        content += "SIM_NSNAP 2 \n"
        content += "SIM_SEED %d \n" % (obsId % 10000 + 4)
        content += "Opsim_rawseeing -1 \n"

        # Add the sky information
        content += "rightascension %.6f \n" % ra
        content += "declination %.6f \n" % dec
        content += "rotskypos %.6f \n" % rot

        if (filePath is not None):
            self.writeToFile(filePath, content=content, mode="w")

        return content

    def getDefaultOpdInstance(self, obsId, aFilterId, filePath=None):
        """
        
        Get the default optical path difference (OPD) instance catalog. It is 
        noted that there is no sky information for OPD.
        
        Arguments:
            obsId {[int]} -- Observation ID.
            aFilterId {[int]} -- PhoSim lsst filter ID.
        
        Keyword Arguments:
            filePath {[str]} -- File path to save the instance. (default: {None})
        
        Returns:
            [str] -- Instance command used in PhoSim.
        """
        
        # Observation Parameters
        content = ""
        content += "Opsim_obshistid %d \n" % obsId
        content += "Opsim_filter %d \n" % aFilterId
        content += "SIM_VISTIME 15.0 \n"
        content += "SIM_NSNAP 1 \n"

        if (filePath is not None):
            self.writeToFile(filePath, content=content, mode="w")

        return content

    def runPhoSim(self, argstring="-v"):
        """
        
        Run the PhoSim program.
        
        Arguments:
            phosimDir {[str]} -- PhoSim directory.
            
        Keyword Arguments:
            argstring {[str]} -- Arguments for PhoSim. (default: {"-h})
        """

        # Path of phosim.py script 
        phosimRunPath = os.path.join(self.phosimDir, "phosim.py")

        # Command to execute the python
        command = " ".join(["python", phosimRunPath])

        # Run the PhoSim with the related arguments
        self.__runProgram(command, argstring=argstring)

    def getPhoSimArgs(self, instance, extraCommand=None, numProc=1, numThread=1, outputDir=None, 
                      instrument="lsst", e2ADC=1, logFilePath=None):
        """
        
        Get the arguments needed to run the PhoSim.
        
        Arguments:
            instance {[str]} -- Instance catalog file.
        
        Keyword Arguments:
            extraCommand {[str]} -- Command file to modify the default physics. (default: {None})
            numProc {int} -- Number of processors. (default: {1})
            numThread {int} -- Number of threads. (default: {1})
            outputDir {[str]} -- Output image directory. (default: {None})
            instrument {str} -- Instrument site directory. (default: {"lsst"})
            e2ADC {int} -- Whether to generate amplifier images (1 = true, 0 = false). (default: {1})
            logFilePath {[str]} -- Log file path for PhoSim calculation log. (default: {None})
        
        Returns:
            [str] -- Arguments to run the PhoSim.
        """
        
        # Use the absolute path
        instance = os.path.abspath(instance)
        
        if (extraCommand is not None):
            extraCommand = os.path.abspath(extraCommand)

        if (outputDir is not None):
            outputDir = os.path.abspath(outputDir)

        if (logFilePath is not None):
            logFilePath = os.path.abspath(logFilePath)

        # Prepare the argument list
        argString = "%s -i %s -e %d" % (instance, instrument, e2ADC)

        if (extraCommand is not None):
            argString += " -c %s" % extraCommand

        if (numProc > 1):
            argString += " -p %d" % numProc

        if (numThread > 1):
            argString += " -t %d" % numThread

        if (outputDir is not None):
            argString += " -o %s" % outputDir

        if (logFilePath is not None):
            argString += " > %s 2>&1" % logFilePath

        return argString

    def __runProgram(self, command, binDir=None, argstring=None):
        """
        
        Run the program w/o arguments.
        
        Arguments:
            command {[string]} -- Command of application.
        
        Keyword Arguments:
            binDir {[str]} -- Directory of binary application. (default: {None})
            argstring {[str]} -- Arguments of program. (default: {None})
        
        Raises:
            RuntimeError -- There is the error in running the program.
        """

        # Directory of binary application
        if (binDir is not None):
            command = os.path.join(binDir, command)

        # Arguments for the program
        if (argstring is not None):
            command += (" " + argstring)

        # Call the program w/o arguments
        if (subprocess.call(command, shell=True) != 0):
            raise RuntimeError("Error running: %s" % command)

    def writeToFile(self, filePath, content=None, sourceFile=None, mode="a"):
        """
        
        Write the file based on the content to put or the source file to copy with.
        
        Arguments:
            filePath {[str]} -- File path to write.
        
        Keyword Arguments:
            content {[str]} -- Content to write into the file. (default: {None})
            sourceFile {[str]} -- Source file to write its content into the file. (default: {None})
            mode {[str]} -- Overwrite ("w") or append ("a") the file. (default: {"a"})
        """

        if mode not in ("w", "a"):
            raise ValueError("Mode: %s is not supported." % mode)

        if (content is not None) or (sourceFile is not None):

            # Open the file. If the file path does not exist, the new file will be generated.
            # Use the append instead of 
            fid = open(filePath, mode)

            # Write the content into the file
            if (content is not None):
                fid.write(content)

            # Write the content of source file into the file
            if (sourceFile is not None):
                fSrc = open(sourceFile, "r")
                fid.write(fSrc.read())
                fSrc.close()

            # Close the file
            fid.close()

if __name__ == "__main__":

    # Directory of phosim
    phosimDir = "/Users/Wolf/Documents/bitbucket/phosim_syseng2"

    # Instantiate the phosim communicator
    phosimCom = PhosimCommu(phosimDir)

    # Run PhoSim
    # phosimCom.runPhoSim()

    # DOF perturbation
    # dof = np.random.rand(50)
    # content = phosimCom.doDofPert(dof)
    # print(content)

    # Get the filter type
    # afilter = "g"
    # afilterId = phosimCom.getFilterId(afilter)
    # print(afilterId)

    # Do the surface perturbation
    # zkInMm = np.random.rand(22)
    # content = phosimCom.doSurfPert(3, zkInMm)
    # print(content)

    # Do the surface map perturbation
    surfId = 3
    surfFilePath = "/data/map1.txt"
    relScale = 1
    content = phosimCom.doSurfMapPert(surfId, surfFilePath, relScale)
    print(content)



