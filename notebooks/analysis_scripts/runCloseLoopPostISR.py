import os
import argparse
import numpy as np

from lsst.ts.wep.WepController import WepController
from lsst.ts.wep.Utility import CamType, FilterType, getModulePath, mapFilterRefToG, DefocalType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory

def main(closed_loop_input_dir, wfs_zer_output, detector_list, postage_img_dir,
         save_postage_stamps=False, run_deblender=True):

    """
    closed_loop_input_dir: string
        The AOS close loop `input` directory.

    wfs_zer_output: string
        The output filename for the wfs Zernikes

    detector_list: list of strings
        The list of detectors with donut images available.

    postage_img_dir: string
        The directory for postage stamps.

    save_postage_stamps: bool, default=False
        Set to save the postage stamps created for each donut

    run_deblender: bool, default=True
        Run the WFS Zernike estimation using donuts that require deblending
    """

    wep_calc = WEPCalculationFactory.getCalculator(CamType.ComCam, closed_loop_input_dir)
    wep_calc.wepCntlr.setPostIsrCcdInputs(os.path.join(closed_loop_input_dir, 'rerun/run1'))
    intraObsId = 9006002
    extraObsId = 9006001
    visitList = [intraObsId, extraObsId]

    neighborStarMap = wep_calc._getTargetStar(visitList, defocalState=DefocalType.Intra)

    isrImgMap = wep_calc.wepCntlr.getPostIsrImgMapByPistonDefocal(detector_list, visitList)

    donut_map = wep_calc.wepCntlr.getDonutMap(neighborStarMap, isrImgMap, FilterType.REF,
                                              doDeblending=run_deblender, postageImg=save_postage_stamps,
                                              postageImgDir=postage_img_dir)

    donutMap = wep_calc.wepCntlr.calcWfErr(donut_map, postage_img_dir)

    listOfWfErr = wep_calc._populateListOfSensorWavefrontData(donutMap)

    zerDict = {}
    for wfErrObj in listOfWfErr:
        zerDict[wfErrObj.getSensorId()] = wfErrObj.getAnnularZernikePoly()

    sensorKeys = list(zerDict.keys())
    sensorKeys.sort()

    zerList = []
    for sensorKey in sensorKeys:
        zerList.append(zerDict[sensorKey])
    zerArray = np.array(zerList)

    np.savetxt(wfs_zer_output, zerArray, header='The followings are ZK in um from z4 to z22:')

if __name__ == '__main__':

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop after ISR is already performed.")
    parser.add_argument("--closed_loop_input_dir", type=str,
                        help="The 'input' directory created by AOS close-loop code.")
    parser.add_argument("--wfs_zer_output", type=str,
                        help="The output filename for the wfs Zernikes")
    parser.add_argument("--postage_img_dir", type=str,
                        help="The directory for postage stamps.")
    parser.add_argument("--save_postage_stamps", default=False, action='store_true',
                        help="Set tag to save the postage stamps created for each donut.")
    parser.add_argument("--run_deblender", default=False, action='store_true',
                        help="Set tag to use donuts that require deblending")
    args = parser.parse_args()

    detector_list = ['R:2,2 S:0,0', 'R:2,2 S:0,1', 'R:2,2 S:0,2']

    main(args.closed_loop_input_dir, args.wfs_zer_output, detector_list,
         args.postage_img_dir, save_postage_stamps=args.save_postage_stamps,
         run_deblender=args.run_deblender)