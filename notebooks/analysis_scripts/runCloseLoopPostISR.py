import os
import argparse
import numpy as np

from lsst.utils import getPackageDir
from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.WepController import WepController
from lsst.ts.wep.Utility import CamType, FilterType, getModulePath, mapFilterRefToG, DefocalType
from lsst.ts.wep.ctrlIntf.WEPCalculationFactory import WEPCalculationFactory

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory
from lsst.ts.wep.ctrlIntf.MapSensorNameAndId import MapSensorNameAndId
from lsst.ts.ofc.ctrlIntf.FWHMSensorData import FWHMSensorData


import sys
sys.path.append('../analysis_tools/')
import analysisTools as at

def main(closed_loop_input_dir, wfs_zer_output, save_postage_stamps=True, run_deblender=True, 
         select_sensor = 'comcam', rerun='run1',
         db_filename = 'bsc1.db3', gaia_field_name='', rotAngInDeg = 0, bscDbType='file',
         setting_filename = 'default.yaml', raInDeg=0, decInDeg=0,expWcs=False):

    """
    closed_loop_input_dir: string
        The AOS close loop `input` directory.

    wfs_zer_output: string, eg. 'wfs.gaia.zer' 
        The output filename for the wfs Zernikes

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
        
    rotAngInDeg : float, rotation angle in degrees - setting passed to wep_calc
    
    bscDbType : str, desired database type, passed on to update the setting file (eg. default.yaml)
         default='file', could be 'refCat', 'image' 
    
    setting_filename: str, name of the setting file for ts_wep/policy/ , eg. 'default.yaml'
    
    """
    
    # update the setting for bscDbType,  if required 
    # for that we need to manually point to ts_wep
    # one way to do it is to use lsst.utils.getPackageDir  
    
    path_to_ts_wep = getPackageDir("ts_wep")
    path_to_setting_file = os.path.join(path_to_ts_wep, 'policy',setting_filename)
    settingFile = ParamReader(filePath=path_to_setting_file)
    bscDbTypeInFile = settingFile.getSetting("bscDbType")
    print('The setting file contains bscDbType : %s'%bscDbTypeInFile)

    if bscDbTypeInFile != bscDbType:
        # In the following we update the setting for bscDbType,
        # saving the change in the default.yaml file 
        settingFile.updateSetting("bscDbType", bscDbType)
        settingFile.saveSetting(filePath=path_to_setting_file)

        # check that the change indeed took place 
        settingFile = ParamReader(filePath=path_to_setting_file)
        print('After change: ', settingFile.getSetting("bscDbType"))
    
    bscDataDir = os.path.join(path_to_ts_wep, 'tests/testData')
    if db_filename  in os.listdir(bscDataDir):
        os.remove(os.path.join(bscDataDir,db_filename))
        print('Removed old %s file'%db_filename)
    
    # update expWcs setting ... 
    #expWcs = False # using PhoSim WCS Sol 
    settingFile= ParamReader(filePath=path_to_setting_file)
    settingFile.updateSetting("expWcs", expWcs)
    settingFile.saveSetting(filePath=path_to_setting_file)
    settingFile = ParamReader(filePath=path_to_setting_file)
    print('After change: expWcs : ', settingFile.getSetting("expWcs"))


 
    ########################
    # Initialize wep_calc
    ########################
    print('select_sensor is %s'%select_sensor)
    if select_sensor == 'comcam':
        wep_calc = WEPCalculationFactory.getCalculator(CamType.ComCam, closed_loop_input_dir)
    elif select_sensor  == 'lsstcam':
        wep_calc = WEPCalculationFactory.getCalculator(CamType.LsstCam, closed_loop_input_dir)
    elif select_sensor == 'lsstfamcam':
        wep_calc = WEPCalculationFactory.getCalculator(CamType.LsstFamCam, closed_loop_input_dir)
    
    #if bscDbType == 'file':

    baseOutputDir = closed_loop_input_dir[:-len('/input')]
    print('\nbaseOutputDir = %s'%baseOutputDir)
    # the iteration directory 
    iterCount = 0
    iterDirName = "%s%d" % ("iter", iterCount)

    # Set the path to sky file directory :   iter0/pert
    if bscDbType == 'file':
        outputDirName = "pert"
        outputDir = os.path.join(baseOutputDir, iterDirName, outputDirName)
        skyInfoFileName = 'skyInfo.txt'
        outputSkyInfoFilePath = os.path.join(outputDir, skyInfoFileName)
        print('\nSetting sky file as %s'%outputSkyInfoFilePath)
        wep_calc.setSkyFile(outputSkyInfoFilePath)

    # update name of bscDbFile ... 
    #settingFilePath = os.path.join(path_to_ts_wep, 'policy/default.yaml')
    #settingFile = ParamReader(filePath=settingFilePath)
    settingFile = wep_calc.getSettingFile()
    dbRelativePath = 'tests/testData/%s'%db_filename
    settingFile.updateSetting("defaultBscPath", dbRelativePath)
    
    # update wep_calc pointing based  on the GAIA field name 
    if len(gaia_field_name)>0:
        raInDeg, decInDeg = at.getRaDecFromGaiaField(gaia_field_name)
    
    # if GAIA field name is not provided, we assume raInDeg, decInDeg as boresight 
    print('\nSetting telescope boresight as (%.2f, %.2f)'%(raInDeg,decInDeg))
    wep_calc.setBoresight(raInDeg, decInDeg)

    print('\nSetting telescope rotation angle as %.2f'%rotAngInDeg)
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

    outputPostageDirName = "postage"
    if save_postage_stamps :
        outputPostageDir  = os.path.join(baseOutputDir,iterDirName,
                            outputPostageDirName)
        if (not os.path.exists(outputPostageDir)):
            os.makedirs(outputPostageDir)
    else:
        outputPostageDir = None


    print('\nGetting the donut map')
    print('run_deblender is ', run_deblender )
    donut_map = wep_calc.wepCntlr.getDonutMap(neighborStarMap, isrImgMap, FilterType.REF,
                                              doDeblending=run_deblender, postageImg=save_postage_stamps,
                                              postageImgDir=outputPostageDir, verbose=True)
    print('\nCalculating wavefront error')
    donutMap = wep_calc.wepCntlr.calcWfErr(donut_map, outputPostageDir, verbose=False)

    # if select_sensor is 'lsstcam' : 
    #     sensorNameToIdFileName='sensorNameToIdWfs.yaml'
    # else:
    #     sensorNameToIdFileName='sensorNameToId.yaml'
    # print('Using sensor to ID translation from %s'%sensorNameToIdFileName)

    # I just make a new dictionary, and since I know from above that 
    # error for R:0,0 S:2,2,A is the same as for R:0,0 S:2,2,B,
    # I keep just one of these 
    if select_sensor == 'lsstcam':
        donutMapAbbrev  = {}
        detectors = list(donutMap.keys())

        for detector in detectors:
            #raft,sensor = parseAbbrevDetectorName(abbrevDetectorName(detector))
            #detect = '%s_%s'%(raft,sensor)
            donutMapAbbrev[detector[:-2]] = donutMap[detector]
        donutMap = donutMapAbbrev.copy()

    listOfWfErr = wep_calc._populateListOfSensorWavefrontData(donutMap)

    
    # Record the wfs error with the same order as OPD for the comparison
    # lines below taken from 
    # phosimCmpt.reorderAndSaveWfErrFile(listOfWfErr, sensorNameList,
    #                                      zkFileName=wfsZkFileName)
    zerDict = {}
    for wfErrObj in listOfWfErr:
        zerDict[wfErrObj.getSensorId()] = wfErrObj.getAnnularZernikePoly()

    sensorKeys = list(zerDict.keys())
    sensorKeys.sort()

    zerList = []
    for sensorKey in sensorKeys:
        zerList.append(zerDict[sensorKey])
    zerArray = np.array(zerList)
    
    # this is the directory where both opd.zer.wfs and arrow.zer.wfs get stored /// 
    # as well as the PSSN.txt ...
    outputImgDir  = os.path.join(baseOutputDir, 'iter0','img')
    wfs_output_path = os.path.join(outputImgDir, wfs_zer_output)
    print('\nSaving the list of wavefront errors as %s'%wfs_output_path)
    np.savetxt( wfs_output_path , zerArray, header='The following are ZK in um from z4 to z22:')


    # calculate corrections with OFC calc...
    # OFC is converting wavefront errors into corrections
    # utilized by M1M3 (figure), M2 (position and figure), and Hexapod (position).
    if select_sensor == 'comcam':
        ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
    if (select_sensor == 'lsstcam' ) or (select_sensor == 'lsstfamcam'): # is this true ? 
        ofcCalc = OFCCalculationFactory.getCalculator(InstName.LSST)

    ofcCalc.setFilter(FilterType.REF)
    ofcCalc.setRotAng(rotAngInDeg)
    ofcCalc.setGainByPSSN()

    # below comes from a call to  PhosimCmpt.py that 
    # returns the list of FWHM per sensorId from the PSSN file 
    #listOfFWHMSensorData = phosimCmpt.getListOfFwhmSensorData(
    #                                         opdPssnFileName, sensorNameList)
    pssnFileName = "PSSN.txt"
    
    filePath = os.path.join(outputImgDir, pssnFileName)
    data = np.loadtxt(filePath)
    fwhmData = data[1, :-1]

    sensorNameToIdFileName  = 'sensorNameToId.yaml'
    mapSensorNameAndId = MapSensorNameAndId(sensorNameToIdFileName)
     
    # if comcam , sensorIdList is just comcam sensors:
    if select_sensor  == 'comcam': 
        refSensorNameList = _getComCamSensorNameList()
        sensorIdList = mapSensorNameAndId.mapSensorNameToId(refSensorNameList)

    # if lsstcam / lsstfamcam,   OPD calc for LsstFamCam, but need to choose 31 of 189 sensors 
    if (select_sensor == 'lsstcam' ) or (select_sensor == 'lsstfamcam'):
        refSensorNameList = _getLsstFamCamSensorNameList()
        
        # getting the names of sensors that 
        # fall within opd field Idx...
        path_to_ts_ofc = getPackageDir("ts_ofc")
        mappingFilePath = os.path.join(path_to_ts_ofc , 'policy/lsst', 'sensorNameToFieldIdx.yaml')
        _mappingFile = ParamReader()
        _mappingFile.setFilePath(mappingFilePath)

        fieldIdx = []
        for sensor in refSensorNameList:
            field = _mappingFile.getSetting(sensor)
            fieldIdx.append(int(field))
            
        uniqueFieldIdx = np.unique(fieldIdx)
        uniqueFieldIdxLt31 = uniqueFieldIdx[uniqueFieldIdx<31]

        #This shows all sensors corresponding to each fieldIdx 
        oneSensorPerFieldIdx = []
        for field in uniqueFieldIdxLt31:
            #print(field, np.array(refSensorNameList)[fieldIdx == field])
            oneSensorPerFieldIdx.append(np.array(refSensorNameList)[fieldIdx == field][0])
            

        # I pick one sensor for each fieldIdx 
        # this list is used to make the listOfFWHMSensorData : 
        sensorIdList = mapSensorNameAndId.mapSensorNameToId(oneSensorPerFieldIdx)

    
    # make a list of fwhm sensor data 
    listOfFWHMSensorData = []
    for sensorId, fwhm in zip(sensorIdList, fwhmData):
        fwhmSensorData = FWHMSensorData(sensorId, np.array([fwhm]))
        listOfFWHMSensorData.append(fwhmSensorData)

    ofcCalc.setFWHMSensorDataOfCam(listOfFWHMSensorData)
    ofcCalc.calculateCorrections(listOfWfErr)
    
    # Save the aggregated degree of freedom to a file 
    dofInUm = ofcCalc.getStateAggregated()
    
    # I skip setting PhoSimCmpt dofInUm , since here we just run one iteration
    # Save the DOF in um data to file for the next iteration
    # parts of 
    # phosimCmpt.saveDofInUmFileForNextIter(
    #        dofInUm, dofInUmFileName=dofInUmFileName)
    dofInUmFileName="dofPertInNextIter.mat"
    outputDir= os.path.join(baseOutputDir,'iter0','pert')
    filePath = os.path.join(outputDir, dofInUmFileName)
    header = "The following are the DOF in um:"
    np.savetxt(filePath, np.transpose(dofInUm), header=header)


def _getComCamSensorNameList():
    chips  = ['00','01','02',
              '10','11','12',
              '20','21','22']
    rafts = ['22']
    sensors = []
    for r in rafts:
        for c in chips:
            s = "R%s_S%s"%(r,c)
            sensors.append(s)
    sensorNameList = sensors
    return sensorNameList

def _getLsstFamCamSensorNameList():
    # I assume it includes the ComCam 
    chips  = ['00','01','02',
              '10','11','12',
              '20','21','22']
    rafts = ['14','24','34', 
        '03','13','23','33','43',
        '02','12','22','32','42',
        '01','11','21','31','41',
             '10','20','30']
    sensors = []
    for r in rafts:
        for c in chips:
            s = "R%s_S%s"%(r,c)
            sensors.append(s)
    sensorNameList = sensors
    return sensorNameList



if __name__ == '__main__':

    # Set the parser
    parser = argparse.ArgumentParser(
        description="Run AOS closed-loop after ISR is already performed.")
    parser.add_argument("--closed_loop_input_dir", type=str,
                        help="The 'input' directory created by AOS close-loop code (absolute path)")
    parser.add_argument("--wfs_zer_output", type=str, default='wfs.out.zer',
                        help="The output filename for the wfs Zernikes (saved inside closed_loop_input_dir)")
    # parser.add_argument("--postage_img_dir", type=str,
    #                     help="The directory for postage stamps (absolute path)")
    parser.add_argument("--save_postage_stamps", default=True, action='store_true',
                        help="Set tag to save the postage stamps created for each donut.")
    parser.add_argument("--run_deblender", default=False, action='store_true',
                        help="Set tag to use donuts that require deblending")
    parser.add_argument("--select_sensor", default='comcam', 
                        help='Choose sensor family: comcam, lsstcam, or lsstfamcam')
    parser.add_argument("--rerun", default='run1', help='The postISR rerun dir')
    parser.add_argument("--db_filename", default='bsc1.db', help='The database filename')
    parser.add_argument("--gaia_field_name", default='', help='The name of GAIA field')
    parser.add_argument("--bscDbType", default='file',
                        help='desired database type: refCat, file, or image' )
    parser.add_argument("--setting_filename", default='default.yaml')
    parser.add_argument("--rotAngInDeg", default=0, help='rotation angle in degrees')
    parser.add_argument("--raInDeg", default=0, help='boresight right ascension in degrees')
    parser.add_argument("--decInDeg", default=0, help='boresight declination in degrees')
    parser.add_argument("--expWcs", default=False, help='whether to use PhosimWcsSol (True) or default WcsSol (False)')
    args = parser.parse_args()

    # detector_list = ['R:2,2 S:0,0', 'R:2,2 S:0,1', 'R:2,2 S:0,2']

    main(closed_loop_input_dir=args.closed_loop_input_dir, wfs_zer_output=args.wfs_zer_output, 
    save_postage_stamps=args.save_postage_stamps, 
    run_deblender=args.run_deblender,select_sensor = args.select_sensor, rerun=args.rerun,
    db_filename = args.db_filename, gaia_field_name=args.gaia_field_name,
        setting_filename=args.setting_filename,
         bscDbType = args.bscDbType,rotAngInDeg=args.rotAngInDeg, raInDeg=args.raInDeg,
         decInDeg=args.decInDeg, expWcs =args.expWcs)