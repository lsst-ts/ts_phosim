import os, unittest
import numpy as np
from scipy.interpolate import Rbf

from cwfs.Tool import ZernikeAnnularFit, ZernikeAnnularEval

from wepPhoSim.MirrorSim import MirrorSim
from wepPhoSim.CoTransform import M1CRS2ZCRS

class M1M3Sim(MirrorSim):
    
    def __init__(self, surf=None, mirrorDataDir=None):
        """
        
        Initiate the M1M3Sim object.
        
        Keyword Arguments:
            surf {[ndarray]} -- Mirror surface along z direction. (default: {None})
            mirrorDataDir {[str]} -- Mirror data directory. (default: {None})
        """

        # Outer radius of M1 mirror in m
        R = 4.180

        # Inner radius of M1 mirror in m
        Ri = 2.558

        # Outer radius of M3 mirror in m
        R3 = 2.508

        # Inner radius of M3 mirror in m
        R3i = 0.550

        super(M1M3Sim, self).__init__((Ri, R3i), (R, R3), surf=surf, mirrorDataDir=mirrorDataDir)

    def getActForce(self, actForceFileName="M1M3_1um_156_force.DAT"):
        """

        Get the mirror actuator forces in N.

        Keyword Arguments:
            actForceFileName {str} -- Actuator force file name. (default: {"M1M3_1um_156_force.DAT"})

        Returns:
            [ndarray] -- Actuator forces in N.
        """

        forceInN = self.getMirrorData(actForceFileName)

        return forceInN

    def getPrintthz(self, zAngleInRadian, preCompElevInRadian=0, FEAzenFileName="M1M3_dxdydz_zenith.txt", 
                    FEAhorFileName="M1M3_dxdydz_horizon.txt", gridFileName="M1M3_1um_156_grid.DAT"):
        """
        
        Get the mirror print in m along z direction in specific zenith angle.
        
        Arguments:
            zAngleInRadian {[float]} -- Zenith angle in radian.
        
        Keyword Arguments:
            preCompElevInRadian {float} -- Pre-compensation elevation angle in radian. 
                                           (default: {0})
            FEAzenFileName {str} -- Finite element analysis (FEA) model data file name in zenith angle. 
                                    (default: {"M1M3_dxdydz_zenith.txt"})                
            FEAhorFileName {str} -- Finite element analysis (FEA) model data file name in horizon angle. 
                                    (default: {"M1M3_dxdydz_horizon.txt"})
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
        
        Returns:
            [ndarray] -- Corrected projection in m along z direction.
        """

        # Data needed to determine gravitational print through
        data = self.getMirrorData(FEAzenFileName)
        zdx = data[:, 0]
        zdy = data[:, 1]
        zdz = data[:, 2]

        data = self.getMirrorData(FEAhorFileName)
        hdx = data[:, 0]
        hdy = data[:, 1]
        hdz = data[:, 2]

        # Do the M1M3 gravitational correction.
        # Map the changes of dx, dy, and dz on a plane for certain zenith angle
        printthxInM = zdx * np.cos(zAngleInRadian) + hdx * np.sin(zAngleInRadian)
        printthyInM = zdy * np.cos(zAngleInRadian) + hdy * np.sin(zAngleInRadian)
        printthzInM = zdz * np.cos(zAngleInRadian) + hdz * np.sin(zAngleInRadian)

        # Get the bending mode information
        idx1, idx3, bx, by, bz = self.__getMirCoor(gridFileName=gridFileName)

        # Calcualte the mirror ideal shape
        zRef = self.__idealShape(bx*1000, by*1000, idx1, idx3)/1000

        # Calcualte the mirror ideal shape with the displacement
        zpRef = self.__idealShape((bx + printthxInM)*1000, (by + printthyInM)*1000, idx1, idx3)/1000

        # Convert printthz into surface sag to get the estimated wavefront error
        # Do the zenith angle correction by the linear approximation with the ideal shape
        printthzInM = printthzInM - (zpRef - zRef)

        # Normalize the coordinate
        R = self.RinM[0]
        Ri = self.RiInM[0]
        normX = bx/R
        normY = by/R
        obs = Ri/R

        # Fit the annular Zernike polynomials z0-z2 (piton, x-tilt, y-tilt)
        zc = ZernikeAnnularFit(printthzInM, normX, normY, 3, obs)

        # Do the estimated wavefront error correction for the mirror projection
        printthzInM -= ZernikeAnnularEval(zc, normX, normY, obs)

        return printthzInM

    def getTempCorr(self, M1M3TBulk, M1M3TxGrad, M1M3TyGrad, M1M3TzGrad, M1M3TrGrad, 
                    FEAfileName="M1M3_thermal_FEA.txt", gridFileName="M1M3_1um_156_grid.DAT"):
        """
        
        Get the mirror print correction in um along z direction for certain temperature gradient.
        
        Arguments:
            M1M3TBulk {[float]} -- Bulk temperature in degree C. (+/-2sigma spans +/-0.8C)
            M1M3TxGrad {[float]} -- Temperature gradient along x direction in degree C. 
                                    (+/-2sigma spans 0.4C)
            M1M3TyGrad {[float]} -- Temperature gradient along y direction in degree C. 
                                    (+/-2sigma spans 0.4C)
            M1M3TzGrad {[float]} -- Temperature gradient along z direction in degree C. 
                                    (+/-2sigma spans 0.1C)
            M1M3TrGrad {[float]} -- Temperature gradient along r direction in degree C. 
                                    (+/-2sigma spans 0.1C)
        
        Keyword Arguments:
            FEAfileName {str} -- Finite element analysis (FEA) model data file name. 
                                 (default: {"M1M3_thermal_FEA.txt"})
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
        
        Returns:
            [ndarray] -- Corrected projection in um along z direction.
        """
        
        # Data needed to determine thermal deformation
        data = self.getMirrorData(FEAfileName, skiprows=1)

        # These are the normalized coordinates

        # In the original XLS file (thermal FEA model), max(x)=164.6060 in,
        # while 4.18m = 164.5669 in for real mirror. These two numbers do not match.
        # The data here has normalized the dimension (x, y) already.
        # n.b. these may not have been normalized correctly, b/c max(tx)=1.0
        tx = data[:, 0]
        ty = data[:, 1]

        # Below are in M1M3 coordinate system, and in micron

        # Do the fitting in the normalized coordinate
        bx, by = self.__getMirCoor(gridFileName=gridFileName)[2:4]
        R = self.RinM[0]
        normX = bx/R
        normY = by/R

        # Fit the bulk
        tbdz = self.__fitData(tx, ty, data[:, 2], normX, normY)

        # Fit the x-grad
        txdz = self.__fitData(tx, ty, data[:, 3], normX, normY)

        # Fit the y-grad
        tydz = self.__fitData(tx, ty, data[:, 4], normX, normY)

        # Fit the z-grad
        tzdz = self.__fitData(tx, ty, data[:, 5], normX, normY)

        # Fit the r-gradß
        trdz = self.__fitData(tx, ty, data[:, 6], normX, normY)

        # Get the temprature correction
        tempCorrInUm = M1M3TBulk*tbdz + M1M3TxGrad*txdz + M1M3TyGrad*tydz + M1M3TzGrad*tzdz + \
                       M1M3TrGrad*trdz

        return tempCorrInUm

    def genMirSurfRandErr(self, zAngleInRadian, LUTfileName="M1M3_LUT.txt", 
                            forceZenFileName="M1M3_force_zenith.txt", 
                            forceHorFileName="M1M3_force_horizon.txt", 
                            forceInflFileName="M1M3_influence_256.txt", 
                            M1M3ForceError=0.05, nzActuator=156, seedNum=0):
        """
        
        Generate the mirror surface random error.
        
        Arguments:
            zAngleInRadian {[float]} -- Zenith angle in radian.
        
        Keyword Arguments:
            LUTfileName {[str]} -- LUT file name. (default: {"M1M3_LUT.txt"})
            forceZenFileName {str} -- File name of actuator forces along zenith direction. 
                                      (default: {"M1M3_force_zenith.txt"})
            forceHorFileName {str} -- File name of actuator forces along horizon direction. 
                                      (default: {"M1M3_force_horizon.txt"})
            forceInflFileName {str} -- Influence matrix of actuator forces. 
                                       (default: {"M1M3_influence_256.txt"})
            M1M3ForceError {float} -- Ratio of actuator force error. (default: {0.05})
            nzActuator {int} -- Number of actuator along z direction. (default: {156})
            seedNum {int} -- Random seed number. (default: {0})
        
        Returns:
            [ndarray] -- Generated mirror surface random error in m.
        """

        # Get the actuator forces in N of M1M3 based on the look-up table (LUT)
        zangleInDeg = zAngleInRadian/np.pi*180
        LUTforce = self.getLUTforce(zangleInDeg, LUTfileName)

        # Add 5% force error (self.M1M3ForceError). This is for iteration 0 only.
        # This means from -5% to +5% of original actuator's force.
        np.random.seed(seedNum)
        nActuator = len(LUTforce)
        myu = (1 + 2*(np.random.rand(nActuator) - 0.5)*M1M3ForceError)*LUTforce

        # Balance forces along z-axis
        # This statement is intentionally to make the force balance.
        myu[nzActuator-1] = np.sum(LUTforce[:nzActuator]) - np.sum(myu[:nzActuator-1])

        # Balance forces along y-axis
        # This statement is intentionally to make the force balance.
        myu[nActuator-1] = np.sum(LUTforce[nzActuator:]) - np.sum(myu[nzActuator:-1])

        # Get the net force along the z-axis
        zf = self.getMirrorData(forceZenFileName)
        hf = self.getMirrorData(forceHorFileName)
        u0 = zf*np.cos(zAngleInRadian) + hf*np.sin(zAngleInRadian)

        # Calculate the random surface
        G = self.getMirrorData(forceInflFileName)
        randSurfInM = G.dot(myu - u0)

        return randSurfInM

    def getMirrorResInMmInZemax(self, gridFileName="M1M3_1um_156_grid.DAT", numTerms=28, 
                                writeZcInMnToFilePath=None):
        """
        
        Get the residue of surface (mirror print along z-axis) in mm after the fitting with spherical
        Zernike polynomials (zk) under the Zemax coordinate.
        
        Keyword Arguments:
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
            numTerms {int} -- Number of Zernike terms to fit. (default: {28})
            writeZcInMnToFilePath {[str]} -- File path to write the fitted zk in mm. (default: {None})
        
        Returns:
            [ndarray] -- Fitted residue in mm after removing the fitted zk terms in Zemax coordinate.
            [ndarray] -- X position in mm in Zemax coordinate.
            [ndarray] -- Y position in mm in Zemax coordinate.
            [ndarray] -- Fitted zk in mm in Zemax coordinate.
        """
        
        # Get the bending mode information
        bx, by, bz = self.__getMirCoor(gridFileName=gridFileName)[2:5]

        # Transform the M1M3 coordinate to Zemax coordinate
        bxInZemax, byInZemax, surfInZemax = M1CRS2ZCRS(bx, by, self.surf)

        # Get the mirror residue and zk in um
        RinM = self.RinM[0]
        resInUmInZemax, zcInUmInZemax = self._MirrorSim__getMirrorResInNormalizedCoor(surfInZemax,
                                                bxInZemax/RinM, byInZemax/RinM, numTerms)

        # Change the unit to mm
        resInMmInZemax = resInUmInZemax * 1e-3
        bxInMmInZemax = bxInZemax * 1e3
        byInMmInZemax = byInZemax * 1e3
        zcInMmInZemax = zcInUmInZemax * 1e-3

        # Save the file of fitted Zk
        if (writeZcInMnToFilePath is not None):
            np.savetxt(writeZcInMnToFilePath, zcInMmInZemax)

        return resInMmInZemax, bxInMmInZemax, byInMmInZemax, zcInMmInZemax

    def writeMirZkAndGridResInZemax(self, resFile=[], surfaceGridN=200, gridFileName="M1M3_1um_156_grid.DAT",
                                    numTerms=28, writeZcInMnToFilePath=None):
        """
        
        Write the grid residue in mm of mirror surface after the fitting with Zk under the Zemax
        coordinate.
        
        Keyword Arguments:
            resFile {[list]} -- File path to save the grid surface residue map ([M1filePath, M3filePath]). 
                                (default: {[]]})
            surfaceGridN {int} -- Surface grid number. (default: {200})
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
            numTerms {int} -- Number of Zernike terms to fit. (default: {28})
            writeZcInMnToFilePath {[str]} -- File path to write the fitted zk in mm. (default: {None})
        
        Returns:
            [str] -- Grid residue map related data of M1.
            [str] -- Grid residue map related data of M3.
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = self.getMirrorResInMmInZemax(gridFileName=gridFileName,
                                                 numTerms=numTerms, writeZcInMnToFilePath=writeZcInMnToFilePath)[0:3]

        # Get the mirror node
        idx1, idx3 = self.__getMirCoor(gridFileName=gridFileName)[0:2]

        # Grid sample map for M1

        # Change the unit from m to mm
        innerRinMm = self.RiInM[0] * 1e3
        outerRinMm = self.RinM[0] * 1e3

        # Get the residue map used in Zemax
        # Content header: (NUM_X_PIXELS, NUM_Y_PIXELS, delta x, delta y)
        # Content: (z, dx, dy, dxdy)
        contentM1 = self._MirrorSim__gridSampInMnInZemax(resInMmInZemax[idx1], bxInMmInZemax[idx1], 
                                                         byInMmInZemax[idx1], innerRinMm, outerRinMm, 
                                                         surfaceGridN, surfaceGridN, resFile=resFile[0])

        # Grid sample map for M3

        # Change the unit from m to mm
        innerRinMm = self.RiInM[1] * 1e3
        outerRinMm = self.RinM[1] * 1e3

        # Get the residue map used in Zemax
        # Content header: (NUM_X_PIXELS, NUM_Y_PIXELS, delta x, delta y)
        # Content: (z, dx, dy, dxdy)
        contentM3 = self._MirrorSim__gridSampInMnInZemax(resInMmInZemax[idx3], bxInMmInZemax[idx3], 
                                                         byInMmInZemax[idx3], innerRinMm, outerRinMm, 
                                                         surfaceGridN, surfaceGridN, resFile=resFile[1])

        return contentM1, contentM3

    def showMirResMap(self, gridFileName="M1M3_1um_156_grid.DAT", numTerms=28, resFile=[], writeToResMapFilePath=[]):
        """
        
        Show the mirror residue map.
        
        Keyword Arguments:
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"})
            numTerms {int} -- Number of Zernike terms to fit. (default: {28})
            resFile {list} -- File path of the grid surface residue map. (default: {[]})
            writeToResMapFilePath {list} -- File path to save the residue map. (default: {[]})
        """

        # Get the residure map
        resInMmInZemax, bxInMmInZemax, byInMmInZemax = self.getMirrorResInMmInZemax(gridFileName=gridFileName,
                                                                                     numTerms=numTerms)[0:3]

        # Get the mirror node
        idx1, idx3 = self.__getMirCoor(gridFileName=gridFileName)[0:2]

        # Show M1 map
        outerRinMm = self.RinM[0] * 1e3
        self._MirrorSim__showResMap(resInMmInZemax[idx1], bxInMmInZemax[idx1], byInMmInZemax[idx1], outerRinMm,
                                    resFile=resFile[0], writeToResMapFilePath=writeToResMapFilePath[0])

        # Show M3 map
        outerRinMm = self.RinM[1] * 1e3
        self._MirrorSim__showResMap(resInMmInZemax[idx3], bxInMmInZemax[idx3], byInMmInZemax[idx3], outerRinMm,
                                    resFile=resFile[1], writeToResMapFilePath=writeToResMapFilePath[1])

    def __getMirCoor(self, gridFileName="M1M3_1um_156_grid.DAT"):
        """
        
        Get the mirror coordinate and node.
        
        Keyword Arguments:
            gridFileName {str} -- File name of bending mode data. (default: {"M1M3_1um_156_grid.DAT"}) 

        Returns:
            [ndarray] -- M1 node.
            [ndarray] -- M3 node.
            [ndarray] -- x coordinate.
            [ndarray] -- y coordinate.
            [ndarray] -- z coordinate.
        """
        
        # Get the bending mode information
        data = self.getMirrorData(gridFileName)
        
        nodeID = data[:, 0].astype("int")
        nodeM1 = (nodeID == 1)
        nodeM3 = (nodeID == 3)
        
        bx = data[:, 1]
        by = data[:, 2]
        bz = data[:, 3:]

        return nodeM1, nodeM3, bx, by, bz

    def __fitData(self, dataX, dataY, data, x, y):
        """

        Fit the data by radial basis function.

        Arguments:
            dataX {[ndarray]} -- Data x.
            dataY {[ndarray]} -- Data y.
            data {[ndarray]} -- Data to fit.
            x {[ndarray]} -- x coordinate.
            y {[ndarray]} -- y coordinate.

        Returns:
            [ndarray] -- Fitted data.
        """

        # Construct the fitting model
        rbfi = Rbf(dataX, dataY, data)

        # Return the fitted data
        return rbfi(x, y)

    def __idealShape(self, xInMm, yInMm, idxM1, idxM3, dr1=0, dr3=0, dk1=0, dk3=0):
        """

        Calculate the ideal shape of mirror along z direction, which is described by a series of
        cylindrically-symmetric aspheric surfaces.

        Arguments:
            xInMm {[ndarray]} -- coordinate x in 1D array in mm.
            yInMm {[ndarray]} -- coordinate y in 1D array in mm.
            idxM1 {[ndarray]} -- M1 node.
            idxM3 {{ndarray}} -- M3 node.

        Keyword Arguments:
            dr1 {float} -- Displacement of r in mirror 1. (default: {0})
            dr3 {float} -- Displacement of r in mirror 3. (default: {0})
            dk1 {float} -- Displacement of kappa (k) in mirror 1. (default: {0})
            dk3 {float} -- Displacement of kappa (k) in mirror 3. (default: {0})

        Returns:
            [float] -- Ideal mirror surface along z direction.
        """

        # M1 optical design
        r1 = -1.9835e4
        k1 = -1.215
        alpha1 = np.zeros((8, 1))
        alpha1[2] = 1.38e-24

        # M3 optical design
        r3 = -8344.5
        k3 = 0.155
        alpha3 = np.zeros((8, 1))
        alpha3[2] = -4.5e-22
        alpha3[3] = -8.15e-30

        # Get the dimension of input xInMm, yInMm
        nr = xInMm.shape
        mr = yInMm.shape
        if (nr != mr):
            print("In the ideal shape calculation, x is [%d] while y is [%d]." % (nr, mr))
            sys.exit()

        # Calculation the curvature (c) and conic constant (kappa)

        # Mirror 1 (M1)
        c1 = 1/(r1 + dr1)
        k1 = k1 + dk1

        # Mirror 3 (M3)
        c3 = 1/(r3 + dr3)
        k3 = k3 + dk3

        # Construct the curvature, kappa, and alpha matrixes for the ideal shape calculation
        cMat = np.zeros(nr)
        cMat[idxM1] = c1
        cMat[idxM3] = c3

        kMat = np.zeros(nr)
        kMat[idxM1] = k1
        kMat[idxM3] = k3

        alphaMat = np.tile(np.zeros(nr), (8, 1))
        for ii in range(8):
            alphaMat[ii, idxM1] = alpha1[ii]
            alphaMat[ii, idxM3] = alpha3[ii]

        # Calculate the radius
        r2 = xInMm**2 + yInMm**2
        r = np.sqrt(r2)

        # Calculate the ideal surface

        # The optical elements of telescopes can often be described by a series of
        # cylindrically-symmetric aspheric surfaces:
        # z(r) = c * r^2/[ 1 + sqrt( 1-(1+k) * c^2 * r^2 ) ] + sum(ai * r^(2*i)) + sum(Aj * Zj)
        # where i = 1-8, j = 1-N

        z0 = cMat * r2 / (1 + np.sqrt(1 - (1 + kMat) * cMat**2 * r2))
        for ii in range(8):
            z0 += alphaMat[ii, :] * r2**(ii+1)

        # M3 vertex offset from M1 vertex, values from Zemax model
        # M3voffset = (233.8 - 233.8 - 900 - 3910.701 - 1345.500 + 1725.701 + 3530.500 + 900 + 233.800)
        M3voffset = 233.8

        # Add the M3 offset (sum(Aj * Zj), j = 1 - N)
        z0[idxM3] = z0[idxM3] + M3voffset

        # In Zemax, z axis points from M1M3 to M2. the reversed direction (z0>0) is needed.
        # That means the direction of M2 to M1M3.
        return -z0

class M1M3SimTest(unittest.TestCase):
    
    """
    Test functions in M1M3Sim.
    """

    def setUp(self):

        # Directory of M1M3 data
        self.mirrorDataDir = os.path.join("..", "data", "M1M3")

    def testFunc(self):

        # Instantiate the M1M3Sim object
        M1M3 = M1M3Sim()
        self.assertEqual(M1M3.RiInM, (2.558, 0.550))
        self.assertEqual(M1M3.RinM, (4.180, 2.508))

        M1M3.setMirrorDataDir(self.mirrorDataDir)

        forceInN = M1M3.getActForce()
        self.assertEqual(forceInN.shape, (156, 156))

        zAngleInDeg = 27.0912
        zAngleInRadian = zAngleInDeg/180*np.pi
        printthzInM = M1M3.getPrintthz(zAngleInRadian)

        ansFilePath = os.path.join("..", "testData", "testM1M3Func", "M1M3printthz.txt")
        ansPrintthzInM = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(printthzInM-ansPrintthzInM)), 1e-10)

        M1M3TBulk = 0.0902
        M1M3TxGrad = -0.0894
        M1M3TyGrad = -0.1973
        M1M3TzGrad = -0.0316
        M1M3TrGrad = 0.0187
        tempCorrInUm = M1M3.getTempCorr(M1M3TBulk, M1M3TxGrad, M1M3TyGrad, M1M3TzGrad, M1M3TrGrad)

        ansFilePath = os.path.join("..", "testData", "testM1M3Func", "M1M3tempCorr.txt")
        ansTempCorrInUm = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(tempCorrInUm-ansTempCorrInUm)), 1e-10)

        iSim = 6
        randSurfInM = M1M3.genMirSurfRandErr(zAngleInRadian, seedNum=iSim)
        
        ansFilePath = os.path.join("..", "testData", "testM1M3Func", "M1M3surfRand.txt")
        ansRandSurfInM = np.loadtxt(ansFilePath)
        self.assertLess(np.sum(np.abs(randSurfInM-ansRandSurfInM)), 1e-10)

        printthzInUm = printthzInM*1e6
        randSurfInUm = randSurfInM*1e6
        mirrorSurfInUm = printthzInUm + randSurfInUm + tempCorrInUm
        M1M3.setSurfAlongZ(mirrorSurfInUm)

        numTerms = 28
        zcInMmInZemax = M1M3.getMirrorResInMmInZemax(numTerms=numTerms)[3]

        ansFilePath = os.path.join("..", "testData", "testM1M3Func", "sim6_M1M3zlist.txt")
        ansZcInUmInZemax = np.loadtxt(ansFilePath)
        ansZcInMmInZemax = ansZcInUmInZemax*1e-3
        self.assertLess(np.sum(np.abs(zcInMmInZemax[0:numTerms]-ansZcInMmInZemax[0:numTerms])), 1e-9)

        resFile1 = os.path.join("..", "output", "M1res.txt")
        resFile3 = os.path.join("..", "output", "M3res.txt")
        resFile = [resFile1, resFile3]
        M1M3.writeMirZkAndGridResInZemax(resFile=resFile, numTerms=numTerms)
        content1 = np.loadtxt(resFile1)
        content3 = np.loadtxt(resFile3)

        ansFilePath1 = os.path.join("..", "testData", "testM1M3Func", "sim6_M1res.txt")
        ansFilePath3 = os.path.join("..", "testData", "testM1M3Func", "sim6_M3res.txt")
        ansContent1 = np.loadtxt(ansFilePath1)
        ansContent3 = np.loadtxt(ansFilePath3)

        self.assertLess(np.sum(np.abs(content1[0,:]-ansContent1[0,:])), 1e-9)
        self.assertLess(np.sum(np.abs(content1[1:,0]-ansContent1[1:,0])), 1e-9)

        self.assertLess(np.sum(np.abs(content3[0,:]-ansContent3[0,:])), 1e-9)
        self.assertLess(np.sum(np.abs(content3[1:,0]-ansContent3[1:,0])), 1e-9)

        writeToResMapFilePath1 = os.path.join("..", "output", "M1resMap.png")
        writeToResMapFilePath3 = os.path.join("..", "output", "M3resMap.png")
        writeToResMapFilePath = [writeToResMapFilePath1, writeToResMapFilePath3]
        M1M3.showMirResMap(numTerms=numTerms, resFile=resFile, writeToResMapFilePath=writeToResMapFilePath)
        self.assertTrue(os.path.isfile(writeToResMapFilePath1))
        self.assertTrue(os.path.isfile(writeToResMapFilePath3))

        os.remove(resFile1)
        os.remove(resFile3)
        os.remove(writeToResMapFilePath1)
        os.remove(writeToResMapFilePath3)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
    