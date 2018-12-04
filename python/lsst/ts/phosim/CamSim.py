import os
import numpy as np


class CamSim(object):

    # Bound of camera body temperature in degree C.
    MIN_TEMP_IN_DEG_C = 2
    MAX_TEMP_IN_DEG_C = 16

    def __init__(self, camTBinDegC=6.5650, camRotInRad=0, camDataDir=""):
        """Initialization of camera simulator class.

        This class is used to correct the camera distortion.

        Parameters
        ----------
        camTBinDegC : float, optional
            Camera body temperature in degree C. (the default is 6.5650.)
        camRotInRad : float, optional
            Camera rotation angle in radian. (the default is 0.)
        camDataDir : str, optional
            Directory of camera distortion data. (the default is "".)
        """

        self.camTBinDegC = self.setBodyTempInDegC(camTBinDegC)
        self.camRotInRad = self.setRotAngInRad(camRotInRad)
        self.camDataDir = camDataDir

    def setCamDataDir(self, camDataDir):
        """Set the camera distortion data directory.

        Parameters
        ----------
        camDataDir : str
            Camera distortion data directory.
        """

        self.camDataDir = camDataDir

    def getCamDataDir(self):
        """Get the camera distortion data directory.

        Returns
        -------
        str
            Camera distortion data directory.
        """

        return self.camDataDir

    def setRotAngInRad(self, rotAngInRad):
        """Set the camera rotation angle in radian.

        The angle should be in (-pi/2, pi/2).

        Parameters
        ----------
        rotAngInRad : float
            Rotation angle in radian.
        """

        if self._valueInRange(rotAngInRad, -np.pi/2, np.pi/2):
            self.camRotInRad = rotAngInRad

    def setRotAngInDeg(self, rotAngInDeg):
        """Set the camera rotation angle in degree.

        The angle should be in (-90, 90).

        Parameters
        ----------
        rotAngInDeg : float
            Rotation angle in degree.
        """

        if self._valueInRange(rotAngInDeg, -90, 90):
            self.camRotInRad = np.deg2rad(rotAngInDeg)

    def setBodyTempInDegC(self, tempInDegC):
        """Set the camera body temperature in degree C.

        Parameters
        ----------
        tempInDegC : float
            Temperature in degree C.
        """

        if self._valueInRange(tempInDegC, self.MIN_TEMP_IN_DEG_C,
                              self.MAX_TEMP_IN_DEG_C):
            self.camTBinDegC = tempInDegC

    def _valueInRange(self, value, lowerBound, upperBound):
        """Check the value is in the range or not.

        Parameters
        ----------
        value : float
            Set value.
        lowerBound : float
            Lower bound.
        upperBound : float
            Upper bound.

        Returns
        -------
        bool
            Ture if the value is in the range.

        Raises
        ------
        ValueError
            The setting value should be in (lowerBound, upperBound).
        """

        valueInRange = False
        if (lowerBound <= value <= upperBound):
            valueInRange = True
        else:
            raise ValueError("The setting value should be in (%.3f, %.3f)."
                             % (lowerBound, upperBound))
        return valueInRange

    def getCamDistortionInMm(self, zAngleInRad, camDistType):
        """Get the camera distortion correction in mm.

        Parameters
        ----------
        zAngleInRad : float
            Zenith angle in radian.
        camDistType : CamDistType
            Camera distortion type.

        Returns
        -------
        numpy.ndarray
            Camera distortion in mm.
        """

        # Get the distortion data
        distType = camDistType.name
        dataFile = os.path.join(self.camDataDir, (distType + ".txt"))
        data = np.loadtxt(dataFile, skiprows=1)

        # Calculate the gravity and temperature distortions
        distortion = self._calcGravityDist(data, zAngleInRad) + \
            self._calcTempDist(data)

        # The order/ index of Zernike corrections by Andy in file is different
        # from PhoSim use. Reorder the correction here for PhoSim to use.
        zidx = [1, 3, 2, 5, 4, 6, 8, 9, 7, 10, 13, 14, 12, 15, 11, 19,
                18, 20, 17, 21, 16, 25, 24, 26, 23, 27, 22, 28]

        # The index of python begins from 0.
        distortion = distortion[[x - 1 for x in zidx]]

        return distortion

    def _calcGravityDist(self, camDistData, zAngleInRad):
        """Calculate the distortion from gravity.

        Parameters
        ----------
        camDistData : numpy.ndarray
            Camera distortion data.
        zAngleInRad : float
            Zenith angle in radian.

        Returns
        -------
        numpy.ndarray
            Distortion from gravity.
        """

        # Pre-compensated elevation angle in radian.
        pre_elev = 0

        # Pre-compensated camera rotation angle in radian.
        pre_camR = 0

        distortion = self._gravityDistFunc(camDistData, zAngleInRad,
                                           self.camRotInRad) - \
            self._gravityDistFunc(camDistData, pre_elev, pre_camR)

        return distortion

    def _gravityDistFunc(self, camDistData, zenithAngle, camRotAngle):
        """Gravity distortion function.

        Parameters
        ----------
        camDistData : numpy.ndarray
            Camera distortion data.
        zenithAngle : float
            Zenith angle.
        camRotAngle : float
            Camera rotation angle.

        Returns
        -------
        numpy.ndarray
            Gravity distortion function.
        """

        distFun = camDistData[0, 3:]*np.cos(zenithAngle) + \
            (camDistData[1, 3:]*np.cos(camRotAngle) +
             camDistData[2, 3:]*np.sin(camRotAngle)) * \
            np.sin(zenithAngle)

        return distFun

    def _calcTempDist(self, camDistData):
        """Calculate the distortion from temperature.

        Parameters
        ----------
        camDistData : numpy.ndarray
            Camera distortion data.

        Returns
        -------
        numpy.ndarray
            Distortion from temperature.
        """

        # List of data:
        # [ze. angle, camRot angle, temp (C), dx (mm), dy (mm), dz (mm),
        #  Rx (rad), Ry (rad), Rz (rad)]
        startTempRowIdx = 3
        endTempRowIdx = 10

        # Do the temperature correction by the simple temperature
        # interpolation/ extrapolation

        # If the temperature is too low, use the lowest listed temperature
        # to do the correction.
        if (self.camTBinDegC <= camDistData[startTempRowIdx, 2]):
            distortion = camDistData[startTempRowIdx, 3:]

        # If the temperature is too high, use the highest listed temperature
        # to do the correction.
        elif (self.camTBinDegC >= camDistData[endTempRowIdx, 2]):
            distortion = camDistData[endTempRowIdx, 3:]

        # Get the correction value by the linear fitting
        else:

            # Find the temperature boundary indexes
            p2 = (camDistData[startTempRowIdx:, 2] >
                  self.camTBinDegC).argmax() + startTempRowIdx
            p1 = p2-1

            # Calculate the linear weighting
            w1 = (camDistData[p2, 2] - self.camTBinDegC) / \
                 (camDistData[p2, 2] - camDistData[p1, 2])
            w2 = 1-w1
            distortion = w1*camDistData[p1, 3:] + w2*camDistData[p2, 3:]

        # Minus the reference temperature correction. There is the problem
        # here.
        # If the pre_temp_cam is not on the data list, this statement will
        # fail/ get nothing.

        # Pre-compensated camera temperature in degree C.
        pre_temp_cam = 0

        distortion -= camDistData[(camDistData[startTempRowIdx:, 2] ==
                                  pre_temp_cam).argmax() + startTempRowIdx, 3:]

        return distortion


if __name__ == "__main__":
    pass
