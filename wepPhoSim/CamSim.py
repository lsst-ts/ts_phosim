import os, unittest
import numpy as np

class CamSim(object):

    def __init__(self, camTBinDegC=6.5650, camRotInRad=0, camDataDir=None):
        """
        
        Initiate the CamSim object. This object is used to correct the camera distortion.
        
        Keyword Arguments:
            camTBinDegC {float} -- Camera body temperature in degree C. (default: {6.5650})
            camRotInRad {float} -- Camera rotation angle in radian. (default: {0})
            camDataDir {[str]} -- Directory of camera distortion data. (default: {None})
        """

        self.camTBinDegC = camTBinDegC
        self.camRotInRad = camRotInRad
        self.camDataDir = camDataDir

    def setCamDataDir(self, camDataDir):
        """
        
        Set the camera data directory.
        
        Arguments:
            camDataDir {[str]} -- Camera data directory.
        """

        self.camDataDir = camDataDir

    def setRotAngInRad(self, rotAngInRad):
        """
        
        Set the camera rotation angle in radian. The value should be in (-pi/2, pi/2).
        
        Arguments:
            rotAngInRad {[float]} -- Rotation angle in radian.
        """

        self.__checkValueInRange(rotAngInRad, -np.pi/2, np.pi/2)
        self.camRotInRad = rotAngInRad

    def setRotAngInDeg(self, rotAngInDeg):
        """
        
        Set the camera rotation angle in degree. The value should be in (-90, 90).
        
        Arguments:
            rotAngInDeg {[float]} -- Rotation angle in degree.
        """

        self.__checkValueInRange(rotAngInDeg, -90, 90)
        
        # Change the unit from degree to radian
        self.camRotInRad = rotAngInDeg/180.0*np.pi

    def setBodyTempInDegC(self, tempInDegC):
        """
        
        Set the camera body temperature in degree C. The value should be in (2, 16).
        
        Arguments:
            tempInDegC {[float]} -- Temperature in degree C.
        """

        self.__checkValueInRange(tempInDegC, 2, 16)
        self.camTBinDegC = tempInDegC

    def __checkValueInRange(self, value, lowerBound, upperBound):
        """
        
        Check the value is in the range or not.
        
        Arguments:
            value {[float]} -- Set value.
            lowerBound {[float]} -- Lower bound.
            upperBound {[float]} -- Upper bound.
        
        Raises:
            ValueError -- Value is not in the range.
        """

        if (value < lowerBound) or (value > upperBound):
            raise ValueError("The setting value should be in (%.3f, %.3f)." % (lowerBound, upperBound))

    def getCamDistortionInMm(self, zAngleInRad, distType, pre_elev=0, pre_camR=0, pre_temp_cam=0):
        """

        Get the camera distortion correction in mm.

        Arguments:
            zAngleInRad {[float]} -- Zenith angle in radian.
            distType {[str]} -- Distortion type.

        Keyword Arguments:
            pre_elev {[float]} -- Pre-compensated elevation angle in radian. (default: {0})
            pre_camR {[float]} -- Pre-compensated camera rotation angle in radian. (default: {0})
            pre_temp_cam {[float]} -- Pre-compensated camera temperature in degree C. (default: {0})

        Returns:
            [ndarray] -- Distortion in mm.
        """

        # Check the distortion type
        fileList = os.listdir(self.camDataDir)
        distTypeList = [os.path.splitext(aFile)[0] for aFile in fileList if not aFile.startswith(".")]

        if distType not in distTypeList:
            raise ValueError("The distortion of '%s' is not supported. Only support: %s" % (distType, distTypeList))

        # Path to camera distortion parameter file
        dataFile = os.path.join(self.camDataDir, (distType + ".txt"))

        # Read the distortion
        data = np.loadtxt(dataFile, skiprows=1)

        # Calculate the distortion (dx, dy, dz, rx, ry, rz)
        # Check with Bo for the math here
        distFun = lambda zenithAngle, camRotAngle: data[0, 3:]*np.cos(zenithAngle) + \
                    ( data[1, 3:]*np.cos(camRotAngle) + data[2, 3:]*np.sin(camRotAngle) )*np.sin(zenithAngle)
        distortion = distFun(zAngleInRad, self.camRotInRad) - distFun(pre_elev, pre_camR)

        # Do the temperature correction by the simple temperature interpolation/ extrapolation
        # If the temperature is too low, use the lowest listed temperature to do the correction.
        # List of data:
        # [ze. angle, camRot angle, temp (C), dx (mm), dy (mm), dz (mm), Rx (rad), Ry (rad), Rz (rad)]
        startTempRowIdx = 3
        endTempRowIdx = 10
        if (self.camTBinDegC <= data[startTempRowIdx, 2]):
            distortion += data[startTempRowIdx, 3:]

        # If the temperature is too high, use the highest listed temperature to do the correction.
        elif (self.camTBinDegC >= data[endTempRowIdx, 2]):
            distortion += data[endTempRowIdx, 3:]

        # Get the correction value by the linear fitting
        else:

            # Find the temperature boundary indexes
            p2 = (data[startTempRowIdx:, 2] > self.camTBinDegC).argmax() + startTempRowIdx
            p1 = p2-1

            # Calculate the linear weighting
            w1 = (data[p2, 2] - self.camTBinDegC) / (data[p2, 2] - data[p1, 2])
            w2 = 1-w1
            distortion += w1*data[p1, 3:] + w2*data[p2, 3:]

        # Minus the reference temperature correction. There is the problem here.
        # If the pre_temp_cam is not on the data list, this statement will fail/ get nothing.
        distortion -= data[(data[startTempRowIdx:, 2] == pre_temp_cam).argmax() + startTempRowIdx, 3:]

        # The order/ index of Zernike corrections by Andy in file is different from PhoSim use.
        # Reorder the correction here for PhoSim to use.
        if (distType[-3:] == "zer"):
            zidx = [1, 3, 2, 5, 4, 6, 8, 9, 7, 10, 13, 14, 12, 15, 11, 19,
                    18, 20, 17, 21, 16, 25, 24, 26, 23, 27, 22, 28]
            # The index of python begins from 0.
            distortion = distortion[[x - 1 for x in zidx]]

        return distortion

class CamSimTest(unittest.TestCase):
    """
    Test functions in CamSim.
    """

    def setUp(self):

        # Directory of camera data
        self.camDataDir = os.path.join("..", "data", "camera")

    def testFunc(self):
        
        # Instantiate the CamSim
        camSim = CamSim()

        camSim.setCamDataDir(self.camDataDir)
        self.assertEqual(camSim.camDataDir, self.camDataDir)

        rotAngInRad = 1.0
        camSim.setRotAngInRad(rotAngInRad)
        self.assertEqual(camSim.camRotInRad, rotAngInRad)

        rotAngInDeg = 30.0
        camSim.setRotAngInDeg(rotAngInDeg)
        self.assertEqual(camSim.camRotInRad, rotAngInDeg/180*np.pi)

        tempInDegC = 10.0
        camSim.setBodyTempInDegC(tempInDegC)
        self.assertEqual(camSim.camTBinDegC, tempInDegC)

        zenithAngleInDeg = 27.0912
        zAngleInRad = zenithAngleInDeg/180*np.pi
        distType = "L1S1zer"
        camSim.setBodyTempInDegC(6.5650)
        camSim.setRotAngInRad(-1.2323)
        distortionInMn = camSim.getCamDistortionInMm(zAngleInRad, distType)

        dataFilePath = os.path.join("..", "testData", "testOpdFunc", "sim6_iter0_pert.cmd")
        distData = np.loadtxt(dataFilePath, skiprows=88, usecols=(1, 2, 3))
        idx = (distData[:,0] == 3)
        absDiff = np.sum(np.abs(distortionInMn - distData[idx,-1]))
        self.assertTrue(absDiff < 1e-10)

if __name__ == "__main__":

    # Do the unit test
    unittest.main()
