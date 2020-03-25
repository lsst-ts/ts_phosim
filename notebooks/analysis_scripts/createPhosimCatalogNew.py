import os
import argparse
import numpy as np

from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import FilterType

from lsst.ts.ofc.Utility import InstName
from lsst.ts.ofc.ctrlIntf.OFCCalculationFactory import OFCCalculationFactory

from lsst.ts.phosim.SkySim import SkySim
from lsst.ts.phosim.OpdMetrology import OpdMetrology
from lsst.ts.phosim.Utility import getAoclcOutputPath, getConfigDir


class createPhosimCatalog():

    def createPhosimCatalog(self, numStars, starMag, numFields, 
                             outputFilePath,
                             addStarsCcdEdge = True, 
                             addStarsInAmp = False,
                             addStarsAcrossAmps = False,
                             starSep = None, magList = None):

        """
        numStars: number of stars per field

        starSep: Minimum separation between stars

        magList: Star magnitudes in catalog

        outputFilePath: Filename for output catalog
        """
        
        # Survey parameters
        surveySettingFilePath = os.path.join(getConfigDir(),
                                            "surveySettings.yaml")
        surveySettings = ParamReader(filePath=surveySettingFilePath)
        filterType = FilterType.fromString(
            surveySettings.getSetting("filterType"))
        raInDeg = surveySettings.getSetting("raInDeg")
        decInDeg = surveySettings.getSetting("decInDeg")
        rotAngInDeg = surveySettings.getSetting("rotAngInDeg")

        ofcCalc = self._prepareOfcCalc(filterType, rotAngInDeg)
        skySim = SkySim()
        metr = OpdMetrology()
        metr.setDefaultComcamGQ()
        
        # here starSep is in percentage of ra span
        if addStarsInAmp :  
            print('Making single amp PhoSim cat\n')
            skySim = self._addStarsInAmp(skySim, metr, numStars,
                                     starSep, magList, numFields)


        if addStarsCcdEdge:
            print('Adding stars close to CCD edge\n')
            skySim = self._addStarsCcdEdge(skySim, metr, numStars,
                starMag, numFields)
            

        # here starSep is in number of amps separating the stars: 
        # minimum is 1, maximum is 7 
        # unlike _addStarsInField()
        # this ensures that the stars do not get created  on the 
        # amplifier edges 
        elif addStarsAcrossAmps:
            skySim = self._addStarsAcrossAmps(skySim, metr, numStars,
                                     starSep, magList, numFields)

        skySim.exportSkyToFile(outputFilePath)


    def _prepareOfcCalc(self, filterType, rotAngInDeg):

        ofcCalc = OFCCalculationFactory.getCalculator(InstName.COMCAM)
        ofcCalc.setFilter(filterType)
        ofcCalc.setRotAng(rotAngInDeg)
        ofcCalc.setGainByPSSN()

        return ofcCalc


    def _addStarsCcdEdge(self, skySim, opdMetr, numStars=5, 
        starMag = 16, numFields=9,raMinOffsetPx = 50,
        raOffsetDeltaPx = 50,decOffsetPx = 380):

        raCenterDegList, declCenterDegList = opdMetr.getFieldXY()

        # limit the list to however many fields we are simulating 
        raCenterDegList = raCenterDegList[:numFields]
        declCenterDegList = declCenterDegList[:numFields]

        # at this point, instead of adding any offset, select single ccd 
        # and add few stars on the edge, with distance d from the edge

        declCcdSpan = abs(declCenterDegList[1]-declCenterDegList[0])  # 0.2347 
        raCcdSpan = abs(np.unique(raCenterDegList)[1] - np.unique(raCenterDegList)[0]) # 0.2347 

        raLeftEdgeList =raCenterDegList-raCcdSpan/2

        #Add offset from the edge :
        raNumPx = 4096
        decNumPx = 4096 

        raOffsetPxList = np.arange(raMinOffsetPx,
            raMinOffsetPx+raOffsetDeltaPx*numStars,
            raOffsetDeltaPx)
        raOffsetDegList = (raOffsetPxList / raNumPx) * raCcdSpan

        
        decOffsetDeg = (decOffsetPx / decNumPx) * declCcdSpan


        raCatalog = []
        decCatalog = []

        count = 0 # this iteartes over all CCDs to consider...
        for raOffsetDeg in raOffsetDegList:
            for ra,dec in zip(raLeftEdgeList+raOffsetDeg, 
                declCenterDegList-count*decOffsetDeg):     
                if ra< 0:
                    ra += 360.0
                raCatalog.append(ra)
                decCatalog.append(dec)
            count += 1 


        for i in range(len(raCatalog)):
            starId = i 
            skySim.addStarByRaDecInDeg(starId, raCatalog[i],
                           decCatalog[i], starMag)
        return skySim



    def _addStarsInAmp(self, skySim, opdMetr, numStars, starSep,
                       starMagList, numFields=1):
        starId = 0
        raInDegList, declInDegList = opdMetr.getFieldXY()

        # at this point, instead of adding any offset, select single amp per ccd 
        # and add 2 stars, centered on the amp,
        # with separation d 

        # the amps are 8 in declination direction
        # and 2 in ra direction 
        # together it's an array of 2x8 = 16 amps 
        raAmpNum = 2 
        declAmpNum = 8 

        declCcdSpan = abs(declInDegList[1]-declInDegList[0])  # 0.2347 
        declAmpSpan = declCcdSpan / declAmpNum

        raCcdSpan = abs(np.unique(raInDegList)[1] - np.unique(raInDegList)[0]) # 0.2347 
        raAmpSpan = raCcdSpan / raAmpNum 

        print('Span of amps in dec,ra is ', declAmpSpan, raAmpSpan)

        # Shift the list to the center of one amp : to the left and up from
        # the CCD center ... #4 counting clockwise 
        raInDegAmpCenter= raInDegList - raAmpSpan/2
        declInDegAmpCenter = declInDegList - declAmpSpan/2



        # Add stars to the catalog 
        # starSep is taken as %of the amplifier span 
        starSepDeg  = 0.01 * starSep * raAmpSpan 

        # calculate the span of ra space : raMin, raMax are still arrays 
        # of values of 9 CCDs...
        raMin =  raInDegAmpCenter - starSepDeg / 2
        raMax =  raInDegAmpCenter + starSepDeg / 2

        raCatalog = []
        decCatalog = []
        for i in range(numFields):
            raPositions = np.linspace(raMin[i], raMax[i],  numStars) 
            decPositions = np.ones(numStars)*declInDegAmpCenter[i]
            
            for j in range(numStars):
                ra = raPositions[j]
                dec = decPositions[j]
                #print('input ra,dec:', ra,dec)
                if ra< 0:
                    ra += 360.0
                #print('corrected ra,dec: ', ra,dec)
                skySim.addStarByRaDecInDeg(starId, ra,
                           dec, starMagList[j])
                starId += 1 
         

        return skySim



