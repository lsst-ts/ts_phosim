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
import numpy as np

from lsst.ts.wep.Utility import FilterType

from lsst.ts.phosim.telescope.TeleFacade import TeleFacade
from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.Utility import getConfigDir, getPhoSimPath, getAoclcOutputPath


def main(phosimDir):

    # Settings
    outputDir = getAoclcOutputPath()
    outputImgDir = os.path.join(outputDir, "img")
    os.makedirs(outputImgDir, exist_ok=True)

    configDir = getConfigDir()
    cmdSettingFile = os.path.join(configDir, "cmdFile", "opdDefault.cmd")
    instSettingFile = os.path.join(configDir, "instFile", "opdDefault.inst")

    # Declare the opd metrology and add the interested field points
    metr = OpdMetrology()
    metr.addFieldXYbyDeg(0, 0)
    metr.addFieldXYbyDeg(0.2, 0.3)

    # Set the Telescope facade class
    tele = TeleFacade()
    tele.addSubSys(addCam=True, addM1M3=True, addM2=True)
    tele.setPhoSimDir(phosimDir)

    # Set the survey parameters
    obsId = 9006000
    filterType = FilterType.REF
    zAngleInDeg = 27.0912
    rotAngInDeg = np.rad2deg(-1.2323)
    tele.setSurveyParam(
        obsId=obsId,
        filterType=filterType,
        zAngleInDeg=zAngleInDeg,
        rotAngInDeg=rotAngInDeg,
    )

    # Generate the perturbation
    iSim = 6
    pertCmdFilePath = tele.writePertBaseOnConfigFile(
        outputDir, seedNum=iSim, saveResMapFig=True
    )

    # Update the telescope degree of freedom with sepecific camera dx
    dofInUm = np.zeros(50)
    dofInUm[6] = 1000
    tele.accDofInUm(dofInUm)

    # Write the physical command file
    cmdFilePath = tele.writeCmdFile(
        outputDir,
        cmdSettingFile=cmdSettingFile,
        pertFilePath=pertCmdFilePath,
        cmdFileName="opd.cmd",
    )

    # Write the instance file
    instFilePath = tele.writeOpdInstFile(
        outputDir, metr, instSettingFile=instSettingFile, instFileName="opd.inst"
    )

    # Get the argument to run the PhoSim
    logFilePath = os.path.join(outputImgDir, "opdPhoSim.log")
    argString = tele.getPhoSimArgs(
        instFilePath,
        extraCommandFile=cmdFilePath,
        numPro=2,
        outputDir=outputImgDir,
        e2ADC=0,
        logFilePath=logFilePath,
    )

    # Run the PhoSim
    tele.runPhoSim(argString)

    # Analyze the OPD fits images
    opdFitsFile = os.path.join(outputImgDir, "opd_%d_0.fits.gz" % obsId)
    zk = metr.getZkFromOpd(opdFitsFile=opdFitsFile)[0]
    print("Zk of OPD_0 is %s." % zk)

    wavelengthInUm = tele.getRefWaveLength() * 1e-3
    pssn = metr.calcPSSN(wavelengthInUm, opdFitsFile=opdFitsFile)
    print("Calculated PSSN is %.4f." % pssn)


if __name__ == "__main__":

    phosimDir = getPhoSimPath()
    main(phosimDir)
