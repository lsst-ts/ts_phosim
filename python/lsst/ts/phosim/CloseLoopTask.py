#!/usr/bin/env python

# This file is part of ts_phosim.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os
import shutil

from lsst.ts.wep.Utility import CamType, FilterType
from lsst.ts.wep.ParamReader import ParamReader

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.Utility import getConfigDir as getConfigDirOfc
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath
from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.PhosimCmpt import PhosimCmpt
from lsst.ts.phosim.PlotUtil import plotFwhmOfIters


class CloseLoopTask(object):
    """Close loop task to run the simulation with PhoSim."""

    def runOpd(
        self,
        inst,
        filterTypeName,
        rotCamInDeg,
        m1m3ForceError,
        numPro,
        iterNum,
        baseOutputDir,
        doErsDirCont,
    ):
        """Run the simulation of OPD.

        OPD: Optical path difference.

        Parameters
        ----------
        inst : str
            Instrument to use: comcam or lsstfam.
        filterTypeName : str
            Filter type name: ref, u, g, r, i, z, or y.
        rotCamInDeg : float
            The camera rotation angle in degree (-90 to 90).
        m1m3ForceError : float
            Ratio of M1M3 actuator force error between 0 and 1.
        numPro : int
            Number of processor to run PhoSim.
        iterNum : int
            Number of closed-loop iteration.
        baseOutputDir : str
            Base output directory.
        doErsDirCont : bool
            Do the erase of the content of base output directory or not.
        """

        # Check the input arguments
        camType, instName = self.getCamTypeAndInstName(inst)
        filterType = self.getFilterType(filterTypeName)
        baseOutputDir = self.checkAndCreateBaseOutputDir(baseOutputDir)

        if doErsDirCont:
            self.eraseDirectoryContent(baseOutputDir)

        # Run the simulation
        self._runSim(
            camType,
            instName,
            filterType,
            rotCamInDeg,
            m1m3ForceError,
            numPro,
            iterNum,
            baseOutputDir,
        )

    def getCamTypeAndInstName(self, inst):
        """Get the camera type and instrument name.

        Parameters
        ----------
        inst : str
            Instrument to use: comcam or lsstfam.

        Returns
        -------
        camType : enum 'CamType' in lsst.ts.wep.Utility
            Camera type.
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.

        Raises
        ------
        ValueError
            This instrument is not supported.
        """

        if inst == "comcam":
            return CamType.ComCam, InstName.COMCAM
        elif inst == "lsstfam":
            return CamType.LsstFamCam, InstName.LSSTFAM
        else:
            raise ValueError(f"This instrument ({inst}) is not supported.")

    def getFilterType(self, filterTypeName):
        """Get the filter type.

        Parameters
        ----------
        filterTypeName : str
            Filter type name: ref, u, g, r, i, z, or y.

        Returns
        -------
        filterType : enum 'FilterType' in lsst.ts.wep.Utility
            Filter type.

        Raises
        ------
        ValueError
            This filter type is not supported.
        """

        if filterTypeName == "ref":
            return FilterType.REF
        elif filterTypeName == "u":
            return FilterType.U
        elif filterTypeName == "g":
            return FilterType.G
        elif filterTypeName == "r":
            return FilterType.R
        elif filterTypeName == "i":
            return FilterType.I
        elif filterTypeName == "z":
            return FilterType.Z
        elif filterTypeName == "y":
            return FilterType.Y
        else:
            raise ValueError(f"This filter type ({filterTypeName}) is not supported.")

    def checkAndCreateBaseOutputDir(self, baseOutputDir):
        """Check and create the base output directory.

        This function will create the directory if it does not exist.

        Parameters
        ----------
        baseOutputDir : str
            Base output directory.

        Returns
        -------
        str
            Base output directory.
        """

        if baseOutputDir == "":
            outputDir = getAoclcOutputPath()
        else:
            outputDir = baseOutputDir
        os.makedirs(outputDir, exist_ok=True)

        return outputDir

    def eraseDirectoryContent(self, targetDir):
        """Erase the directory content.

        Parameters
        ----------
        targetDir : str
            Target directory.
        """

        for theFile in os.listdir(targetDir):
            filePath = os.path.join(targetDir, theFile)
            if os.path.isfile(filePath):
                os.unlink(filePath)
            elif os.path.isdir(filePath):
                shutil.rmtree(filePath)

    def _runSim(
        self,
        camType,
        instName,
        filterType,
        rotCamInDeg,
        m1m3ForceError,
        numPro,
        iterNum,
        baseOutputDir,
    ):
        """Run the simulation.

        Parameters
        ----------
        camType : enum 'CamType' in lsst.ts.wep.Utility
            Camera type.
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        filterType : enum 'FilterType' in lsst.ts.wep.Utility
            Filter type.
        rotCamInDeg : float
            The camera rotation angle in degree (-90 to 90).
        m1m3ForceError : float
            Ratio of M1M3 actuator force error between 0 and 1.
        numPro : int
            Number of processor to run PhoSim.
        iterNum : int
            Number of closed-loop iteration.
        baseOutputDir : str
            Base output directory.
        """

        # Prepare the components
        phosimDir = getPhoSimPath()
        phosimCmpt = self._preparePhosimCmpt(
            phosimDir, filterType, rotCamInDeg, m1m3ForceError, numPro
        )

        ofcCalc = self._prepareOfcCalc(instName, filterType, rotCamInDeg)

        # Set the telescope state to be the same as the OFC
        state0 = ofcCalc.getStateAggregated()
        phosimCmpt.setDofInUm(state0)

        # Get the list of referenced sensor name (field positions)
        refSensorNameList = self.getSensorNameListOfFields(instName)

        # Do the iteration
        obsId = 9006000
        opdZkFileName = "opd.zer"
        opdPssnFileName = "PSSN.txt"
        outputDirName = "pert"
        outputImgDirName = "img"
        iterDefaultDirName = "iter"
        dofInUmFileName = "dofPertInNextIter.mat"
        for iterCount in range(iterNum):

            # Set the observation Id
            phosimCmpt.setSurveyParam(obsId=obsId)

            # The iteration directory
            iterDirName = "%s%d" % (iterDefaultDirName, iterCount)

            # Set the output directory
            outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
            phosimCmpt.setOutputDir(outputDir)

            # Set the output image directory
            outputImgDir = os.path.join(baseOutputDir, iterDirName, outputImgDirName)
            phosimCmpt.setOutputImgDir(outputImgDir)

            # Generate the OPD image
            argString = phosimCmpt.getOpdArgsAndFilesForPhoSim(instName)
            phosimCmpt.runPhoSim(argString)

            # Analyze the OPD data
            # Rotate OPD in the reversed direction of camera
            phosimCmpt.analyzeOpdData(
                instName,
                zkFileName=opdZkFileName,
                rotOpdInDeg=-rotCamInDeg,
                pssnFileName=opdPssnFileName,
            )

            # Get the PSSN from file
            pssn = phosimCmpt.getOpdPssnFromFile(opdPssnFileName)
            print("Calculated PSSN is %s." % pssn)

            # Get the GQ effective FWHM from file
            gqEffFwhm = phosimCmpt.getOpdGqEffFwhmFromFile(opdPssnFileName)
            print("GQ effective FWHM is %.4f." % gqEffFwhm)

            # Set the FWHM data
            listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
                opdPssnFileName, refSensorNameList
            )
            ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)

            # Simulate to get the wavefront sensor data from WEP and calculate
            # the DOF
            listOfWfErr = phosimCmpt.mapOpdDataToListOfWfErr(
                opdZkFileName, refSensorNameList
            )
            ofcCalc.calculateCorrections(listOfWfErr)

            # Set the new aggregated DOF to phosimCmpt
            dofInUm = ofcCalc.getStateAggregated()
            phosimCmpt.setDofInUm(dofInUm)

            # Save the DOF file
            phosimCmpt.saveDofInUmFileForNextIter(
                dofInUm, dofInUmFileName=dofInUmFileName
            )

            # Add the observation ID by 1
            obsId += 1

        # Summarize the FWHM
        pssnFiles = [
            os.path.join(
                baseOutputDir,
                "%s%d" % (iterDefaultDirName, num),
                outputImgDirName,
                opdPssnFileName,
            )
            for num in range(iterNum)
        ]
        saveToFilePath = os.path.join(baseOutputDir, "fwhmIters.png")
        plotFwhmOfIters(pssnFiles, saveToFilePath=saveToFilePath)

    def _preparePhosimCmpt(
        self, phosimDir, filterType, rotAngInDeg, m1m3ForceError, numPro
    ):
        """Prepare the component to use PhoSim.

        Parameters
        ----------
        phosimDir : str
            PhoSim path.
        filterType : enum 'FilterType' in lsst.ts.wep.Utility
            Filter type.
        rotAngInDeg : float
            The camera rotation angle in degree (-90 to 90).
        m1m3ForceError : float
            Ratio of M1M3 actuator force error between 0 and 1.
        numPro : int
            Number of processor to run PhoSim.

        Returns
        -------
        PhosimCmpt
            PhoSim component.
        """

        # Set the Telescope facade class
        tele = TeleFacade()
        tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
        tele.setPhoSimDir(phosimDir)

        # Prepare the phosim component
        phosimCmpt = PhosimCmpt(tele)

        # Set the telescope survey parameters
        boresight = (0, 0)
        zAngleInDeg = 27.0912
        phosimCmpt.setSurveyParam(
            filterType=filterType,
            boresight=boresight,
            zAngleInDeg=zAngleInDeg,
            rotAngInDeg=rotAngInDeg,
        )

        # Set the M1M3 force error
        phosimCmpt.setM1M3ForceError(m1m3ForceError)

        # Update the number of processor if necessary
        if numPro > 1:
            settingFile = phosimCmpt.getSettingFile()
            settingFile.updateSetting("numPro", numPro)

        # Set the seed number for M1M3 surface
        seedNum = 6
        phosimCmpt.setSeedNum(seedNum)

        return phosimCmpt

    def _prepareOfcCalc(self, instName, filterType, rotAngInDeg):
        """Prepare the OFC calculator.

        OFC: Optical feedback calculator.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.
        filterType : enum 'FilterType' in lsst.ts.wep.Utility
            Filter type.
        rotAngInDeg : float
            The camera rotation angle in degree (-90 to 90).

        Returns
        -------
        OFCCalculation child in lsst.ts.ofc.ctrlIntf
            Concrete child class of OFCCalculation class.
        """

        ofcCalc = OFCCalculationFactory.getCalculator(instName)

        ofcCalc.setFilter(filterType)
        ofcCalc.setRotAng(rotAngInDeg)
        ofcCalc.setGainByPSSN()

        return ofcCalc

    def getSensorNameListOfFields(self, instName):
        """Get the list of sensor name of fields.

        The list will be sorted based on the field index.

        Parameters
        ----------
        instName : enum 'InstName' in lsst.ts.ofc.Utility
            Instrument name.

        Returns
        -------
        list[str]
            List of sensor name.

        Raises
        ------
        ValueError
            This instrument name is not supported.
        """

        # Check the input
        if instName not in (InstName.COMCAM, InstName.LSSTFAM):
            raise ValueError(f"This instrument name ({instName}) is not supported.")

        # Get the data of sensor name and field index
        filePath = os.path.join(
            getConfigDirOfc(), instName.name.lower(), "sensorNameToFieldIdx.yaml"
        )
        paramReader = ParamReader(filePath=filePath)
        data = paramReader.getContent()

        # Sort the list of sensor name based on field index from low to high
        sensorNameList = [
            senserName
            for senserName, _ in sorted(data.items(), key=lambda item: item[1])
        ]

        return sensorNameList

    @staticmethod
    def setDefaultParser(parser):
        """Set the default parser.

        Parameters
        ----------
        parser : argparse.ArgumentParser
            Input parser.

        Returns
        -------
        argparse.ArgumentParser
            Updated parser.
        """

        parser.add_argument(
            "--inst",
            type=str,
            default="comcam",
            help="Instrument to use: comcam or lsstfam. (default: comcam)",
        )

        parser.add_argument(
            "--filterType",
            type=str,
            default="ref",
            help="Filter type to use: ref, u, g, r, i, z, or y. (default: ref)",
        )

        parser.add_argument(
            "--rotCam",
            type=float,
            default=0.0,
            help="Rotate camera (degree) in counter-clockwise direction. (default: 0.0)",
        )

        parser.add_argument(
            "--m1m3FErr",
            type=float,
            default=0.05,
            help="Ratio of M1M3 actuator force error between 0 and 1. (default: 0.05)",
        )

        parser.add_argument(
            "--numOfProc",
            type=int,
            default=1,
            help="Number of processor to run PhoSim. (default: 1)",
        )

        parser.add_argument(
            "--iterNum",
            type=int,
            default=5,
            help="Number of closed-loop iteration. (default: 5)",
        )

        parser.add_argument("--output", type=str, default="", help="Output directory.")

        parser.add_argument(
            "--clobber",
            default=False,
            action="store_true",
            help="Delete existing output directory.",
        )

        # Return the parser intentionally to emphasise this object has been
        # updated.
        return parser
