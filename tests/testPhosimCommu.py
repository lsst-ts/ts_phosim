import os
import numpy as np
import unittest

from lsst.ts.phosim.PhosimCommu import PhosimCommu
from lsst.ts.phosim.Utility import getModulePath


class TestPhosimCommu(unittest.TestCase):
    """ Test the PhosimCommu class."""

    def setUp(self):

        self.phosimCom = PhosimCommu()

    def testSetAndGetPhoSimDir(self):

        phosimDir = "NotPhoSimDir"
        self.phosimCom.setPhoSimDir(phosimDir)
        self.assertEqual(self.phosimCom.getPhoSimDir(), phosimDir)

    def testFunc(self):
        
        phosimFilterID = self.phosimCom.getFilterId("r")
        self.assertEqual(phosimFilterID, 2)

        surfaceID = self.phosimCom.getSurfaceId("M3")
        self.assertEqual(surfaceID, 2)

        dofInUm = np.zeros(50)
        content = self.phosimCom.doDofPert(dofInUm)
        self.assertEqual(len(content.split("\n")), 51)

        surfId = 1
        zkInMm = [1]
        content = self.phosimCom.doSurfPert(surfId, zkInMm)
        ansContent = "izernike 1 0 1 \n"
        self.assertEqual(content, ansContent)

        surfId = 1
        surfFilePath = "temp.txt"
        surfFilePath = os.path.abspath(surfFilePath)
        relScale = 1
        content = self.phosimCom.doSurfMapPert(surfId, surfFilePath, relScale)
        ansContent = "surfacemap 1 %s 1 \n" % surfFilePath
        self.assertEqual(content, ansContent)

        content = self.phosimCom.doCameraConfig(guidSensorOn=True)
        ansContent = "camconfig 4 \n"
        self.assertEqual(content, ansContent)

        linkSurfId1 = 1
        linkSurfId2 = 2
        content = self.phosimCom.doSurfLink(linkSurfId1, linkSurfId2)
        ansContent = "surfacelink 1 2 \n"
        self.assertEqual(content, ansContent)

        opdId = 0
        fieldXInDeg = 1.0
        fieldYInDeg = 2.0
        wavelengthInNm = 500
        content = self.phosimCom.generateOpd(opdId, fieldXInDeg, fieldYInDeg, wavelengthInNm)
        ansContent = "opd  0\t 1.000000\t 2.000000 500.0 \n"
        self.assertEqual(content, ansContent)

        starId = 0
        ra = 1.0
        dec = 1.0
        magNorm = 2.0
        sedName = "flat.txt"
        content = self.phosimCom.generateStar(starId, ra, dec, magNorm, sedName)
        ansContent = "object  0\t 1.000000\t 1.000000  2.000000 ../sky/flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 star 0.0 none none \n"
        self.assertEqual(content, ansContent)

        obsId = 100
        aFilterId = 1
        content = self.phosimCom.getStarInstance(obsId, aFilterId)
        self.assertEqual(len(content.split("\n")), 8)

        content = self.phosimCom.getOpdInstance(obsId, aFilterId)
        self.assertEqual(len(content.split("\n")), 3)

        instFileName = os.path.join(getModulePath(), "tests", "testData", "temp.inst")
        self.phosimCom.writeToFile(instFileName, content="temp", mode="w")
        argString = self.phosimCom.getPhoSimArgs(instFileName, sensorName="R22_S11")

        absFilePath = os.path.abspath(instFileName)
        ansArgString = "%s -i lsst -e 1 -s R22_S11" % absFilePath

        self.assertEqual(argString, ansArgString)
        os.remove(absFilePath)

        wavelengthInNm = 300.0
        try:
            sedFilePath = self.phosimCom.writeSedFile(wavelengthInNm)
            os.remove(sedFilePath)
        except FileNotFoundError:
            print("Do not find the sed file.")

        try:
            self.phosimCom.runPhoSim()
        except RuntimeError:
            print("Do not find PhoSim directory.")


if __name__ == "__main__":

    # Run the unit test
    unittest.main()
