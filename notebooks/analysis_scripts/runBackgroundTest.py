# import sys
# path_to_comcamCloseLoop = '/astro/store/epyc/projects/lsst_comm/ts_phosim/bin.src/'
# sys.path.append(path_to_comcamCloseLoop)
import os
import argparse
import numpy as np
from baseComcamLoop import baseComcamLoop as comcamLoop
from baseComcamLoop import _eraseFolderContent
from createPhosimCatalog import createPhosimCatalog
from lsst.ts.phosim.Utility import getPhoSimPath, getAoclcOutputPath, getConfigDir


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--numOfProc", type=int, default=1,
                        help="number of processor to run PhoSim (default: 1)")
    parser.add_argument("--iterNum", type=int, default=5,
                        help="number of closed-loop iteration (default: 5)")
    parser.add_argument("--output", type=str, default="",
                        help="output directory")
    parser.add_argument("--testLabel", type=str, default="background")
    parser.add_argument("--testOutput", type=str, default="")
    parser.add_argument("--skyFile", type=str, default="starCat.txt")
    parser.add_argument("--raShift", type=float, default=0.0)
    parser.add_argument("--decShift", type=float, default=0.0)
    parser.add_argument("--opd", default=True, action='store_false')
    parser.add_argument("--defocalImg", default=True, action='store_false')
    parser.add_argument("--flats", default=True, action='store_false')
    parser.add_argument('--phosimCmdFileSuffix', type=str, default="")
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()

    if (args.output == ""):
        outputDir = getAoclcOutputPath()
    else:
        outputDir = args.output
        os.makedirs(outputDir, exist_ok=True)

    testLabel = args.testLabel
    skyFilePath = args.skyFile

    if (args.testOutput == ""):
        testOutputDir = os.path.dirname(os.path.realpath(__file__))
    else:
        testOutputDir = args.testOutput

    os.environ["closeLoopTestDir"] = testOutputDir

    if args.phosimCmdFileSuffix == "":
        opdCmdSettingsFile = 'opdDefault.cmd'
        comcamCmdSettingsFile = 'starDefault.cmd'
    else:
        opdCmdSettingsFile = 'opd%s.cmd' % args.phosimCmdFileSuffix
        comcamCmdSettingsFile = 'star%s.cmd' % args.phosimCmdFileSuffix

    # Clobber
    if args.opd is True:
        _eraseFolderContent(outputDir)
    else:
        if args.flats is True:
            _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
            _eraseFolderContent(os.path.join(outputDir, 'input'))     
        if args.defocalImg is True:
            _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'intra'))
            _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'extra'))

    createCat = createPhosimCatalog()
    raShift = (args.raShift * .2) / 3600 # Convert to degrees
    decShift = (args.decShift * .2) / 3600 # Convert to degrees
    createCat.createPhosimCatalog(1, 0, [15], raShift, decShift,
                                    skyFilePath, numFields=3)

    ccLoop = comcamLoop()
    ccLoop.main(phosimDir, args.numOfProc, args.iterNum, outputDir, args.testLabel,
                isEimg=False, genOpd=args.opd, genDefocalImg=args.defocalImg, 
                genFlats=args.flats, useMinDofIdx=False,
                inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
                opdCmdSettingsFile=opdCmdSettingsFile,
                comcamCmdSettingsFile=comcamCmdSettingsFile)
