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
from lsst.ts.wep.Utility import runProgram


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--numOfProc", type=int, default=1,
                        help="number of processor to run PhoSim (default: 1)")
    parser.add_argument("--testLabel", type=str, default="mag")
    parser.add_argument("--testOutput", type=str, default="")
    parser.add_argument("--skyFile", type=str, default="starCat.txt")
    parser.add_argument("--raShift", type=float, default=0.0)
    parser.add_argument("--decShift", type=float, default=0.0)
    # parser.add_argument("--opd", default=True, action='store_false')
    # parser.add_argument("--defocalImg", default=True, action='store_false')
    # parser.add_argument("--flats", default=True, action='store_false')
    parser.add_argument('--phosimCmdFileSuffix', type=str, default="")
    args = parser.parse_args()

    # Load directory paths
    phosimDir = getPhoSimPath()
    outputDir = getAoclcOutputPath()
    testLabel = args.testLabel
    skyFilePath = args.skyFile
    genFlats = True
    genOpd = True
    genDefocalImg = True

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

    save_images = True
    first_mag = None

    for magVal in np.arange(16.0, 11.1, -0.5):

        if save_images is True:
            dir_name = 'test_output/mag_quickbg_%.2f_output' % magVal
            if os.path.exists(dir_name):
                print('Output Directory exists')
                _eraseFolderContent(dir_name)
            else:
                os.makedirs(dir_name)
            outputDir = dir_name

        # # Clobber
        # if args.opd is True:
        if save_images is False:
            _eraseFolderContent(outputDir)
        # else:
        #     if args.flats is True:
        #         _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
        #         _eraseFolderContent(os.path.join(outputDir, 'input'))     
        #     if args.defocalImg is True:
        #         _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'intra'))
        #         _eraseFolderContent(os.path.join(outputDir, 'iter0', 'img', 'extra'))

        if not genOpd and not genFlats: 
            print('Copying content of /mag_%s/ to re-use the flats and OPD files...' % first_mag)
            argString = '-a test_output/mag_quickbg_%s_output/. %s/' % (first_mag, outputDir)
            runProgram("cp", argstring=argString)

            # ensure that input/raw and input/rerun are empty 
            print('Deleting content of input/raw/ and input/rerun/')
            _eraseFolderContent(os.path.join(outputDir, 'input','raw'))
            _eraseFolderContent(os.path.join(outputDir, 'input','rerun'))

            # remove files that are remade
            argString = os.path.join(outputDir, 'input')
            runProgram("rm", argstring=argString+'/isr*')
            runProgram("rm", argstring=argString+'/registry*')
            runProgram("rm", argstring=argString+'/_mappe*')

        # Clobber
        if genOpd is True:
            print('We will make new OPD files in this run')
            _eraseFolderContent(outputDir)
        else:
            if genFlats is True:
                print('We will make new flats in this run')
                _eraseFolderContent(os.path.join(outputDir, 'fake_flats'))
                _eraseFolderContent(os.path.join(outputDir, 'input'))     
            if genDefocalImg is True:
                print('We will make new defocal images in this run ')
                intraPath = os.path.join(outputDir, 'iter0', 'img', 'intra')
                extraPath = os.path.join(outputDir, 'iter0', 'img', 'extra')
                if os.path.exists(intraPath):
                    _eraseFolderContent(intraPath)
                if os.path.exists(extraPath):
                    _eraseFolderContent(extraPath)

        if first_mag is None:
            first_mag = '%.2f' % magVal

        createCat = createPhosimCatalog()
        raShift = (args.raShift * .2) / 3600 # Convert to degrees
        decShift = (args.decShift * .2) / 3600 # Convert to degrees
        createCat.createPhosimCatalog(1, 0, [magVal], raShift, decShift,
                                      skyFilePath, numFields=3)

        ccLoop = comcamLoop()
        ccLoop.main(phosimDir, args.numOfProc, 1, 
                    outputDir, '%s.%.1f' % (testLabel, magVal), 
                    isEimg=False, genOpd=genOpd, genDefocalImg=True, 
                    genFlats=genFlats, useMinDofIdx=False,
                    inputSkyFilePath=skyFilePath, m1m3ForceError=0.05,
                    opdCmdSettingsFile=opdCmdSettingsFile,
                    comcamCmdSettingsFile=comcamCmdSettingsFile)

        # # Once the necessary data is created we don't need to recreate on every iteration
        genFlats = False
        genOpd = False
        # args.defocalImg = False
