import os
import argparse
import numpy as np

from lsst.ts.wep.WepController import WepController
from lsst.ts.wep.Utility import CamType, FilterType, getModulePath, mapFilterRefToG, DefocalType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory

import sys
sys.path.append('../analysis_tools/')
import analysisTools as at

def main(closed_loop_input_dir, wfs_zer_output, postage_img_dir,
         save_postage_stamps=True, run_deblender=True,select_sensor = 'comcam', rerun='run1',
         db_filename = 'bsc1.db3', gaia_field_name='med', rotAngInDeg = 0,
          ):

    """
    closed_loop_input_dir: string
        The AOS close loop `input` directory.

    wfs_zer_output: string, eg. 'wfs.gaia.zer' 
        The output filename for the wfs Zernikes

    postage_img_dir: string
        The directory for postage stamps.

    save_postage_stamps: bool, default=False
        Set to save the postage stamps created for each donut

    run_deblender: bool, default=True
        Run the WFS Zernike estimation using donuts that require deblending

    select_sensor: str, default='comcam', run on ComCam data, 'lsstcam' - corner wfs sensor data,
        or 'lsstfamcam' - the full array mode data 

    rerun : str, default = 'run1' - which ISR rerun to use ? Useful if have 
        multiple reruns in the same repo 
    
    db_filename: str, default = 'bsc1.db3', may want to use a different one if eg. a previous loop didn't finish,
        or there is a concurrent loop running that may use the same file, to avoid conflicts

    gaia_field_name: str, default = '', could be 'low', 'high', 'Baade', 'Pleiades' - corresponds 
        to the field name of GAIA pointing 
    """
    print('select_sensor is %s'%select_sensor)
    if select_sensor == 'comcam':
        wep_calc = WEPCalculationFactory.getCalculator(CamType.ComCam, closed_loop_input_dir)
    elif select_sensor  == 'lsstcam':
        wep_calc = WEPCalculationFactory.getCalculator(CamType.LsstCam, closed_loop_input_dir)
    elif select_sensor == 'lsstfamcam':
        wep_calc = WEPCalculationFactory.getCalculator(CamType.LsstFamCam, closed_loop_input_dir)
    
    # update name of bscDbFile ... 
    #settingFilePath = os.path.join(path_to_ts_wep, 'policy/default.yaml')
    #settingFile = ParamReader(filePath=settingFilePath)
    settingFile = wep_calc.getSettingFile()
    dbRelativePath = 'tests/testData/%s'%db_filename
    settingFile.updateSetting("defaultBscPath", dbRelativePath)
    
    # update wep_calc pointing 
    if len(gaia_field_name)>0:
        raInDeg, decInDeg = at.getRaDecFromGaiaField(gaia_field_name)
        wep_calc.setBoresight(raInDeg, decInDeg)
    wep_calc.setRotAng(rotAngInDeg)

    wep_calc.wepCntlr.setPostIsrCcdInputs(os.path.join(closed_loop_input_dir, 'rerun/%s'%rerun))
    intraObsId = 9006002
    extraObsId = 9006001
    if select_sensor in ['lsstfamcam','comcam']:
        obsIdList = [intraObsId, extraObsId]
    elif (select_sensor == 'lsstcam'):
        obsIdList = [intraObsId]

    print('\nGetting  neighborStarMap from target stars')
    neighborStarMap = wep_calc._getTargetStar(obsIdList, defocalState=DefocalType.Intra)
    detector_list = list(neighborStarMap)
    
    print('\nGetting postISR image map')
    if len(obsIdList) == 2 : # ComCam or LsstFamCam case
        isrImgMap = wep_calc.wepCntlr.getPostIsrImgMapByPistonDefocal(detector_list, obsIdList)

    elif len(obsIdList) == 1 : #corner WFS case
        isrImgMap = wep_calc.wepCntlr.getPostIsrImgMapOnCornerWfs(detector_list, obsIdList[0])

    print('\nGetting the donut map')
    donut_map = wep_calc.wepCntlr.getDonutMap(neighborStarMap, isrImgMap, FilterType.REF,
                                              doDeblending=run_deblender, postageImg=save_postage_stamps,
                                              postageImgDir=postage_img_dir)
    print('\nCalculating wavefront error')
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
    
    wfs_output_path = os.path.join(closed_loop_input_dir, wfs_zer_output)
    print('\nSaving the list of wavefront errors as %s'%wfs_output_path)
    np.savetxt( wfs_output_path , zerArray, header='The following are ZK in um from z4 to z22:')

if __name__ == '__main__':

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop after ISR is already performed.")
    parser.add_argument("--closed_loop_input_dir", type=str,
                        help="The 'input' directory created by AOS close-loop code (absolute path)")
    parser.add_argument("--wfs_zer_output", type=str,
                        help="The output filename for the wfs Zernikes (saved inside closed_loop_input_dir)")
    parser.add_argument("--postage_img_dir", type=str,
                        help="The directory for postage stamps (absolute path)")
    parser.add_argument("--save_postage_stamps", default=True, action='store_true',
                        help="Set tag to save the postage stamps created for each donut.")
    parser.add_argument("--run_deblender", default=True, action='store_true',
                        help="Set tag to use donuts that require deblending")
    parser.add_argument("--select_sensor", default='comcam', help='Choose sensor family: comcam, lsstcam, or lsstfamcam')
    parser.add_argument("--rerun", default='run1', help='The postISR rerun dir')
    parser.add_argument("--db_filename", default='bsc1.db', help='The database filename')
    parser.add_argument("--gaia_field_name", default='med', help='The name of GAIA field')
    args = parser.parse_args()

    # detector_list = ['R:2,2 S:0,0', 'R:2,2 S:0,1', 'R:2,2 S:0,2']

    main(closed_loop_input_dir=args.closed_loop_input_dir, wfs_zer_output=args.wfs_zer_output, 
    postage_img_dir=args.postage_img_dir, save_postage_stamps=args.save_postage_stamps, 
    run_deblender=args.run_deblender,select_sensor = args.select_sensor, rerun=args.rerun,
    db_filename = args.db_filename, gaia_field_name=args.gaia_field_name)