# setup obs_lsst -t w_2019_44
import pandas as pd
from lsst.obs.base import createInitialSkyWcs
import lsst.obs.lsst as obs_lsst
import lsst.geom as lsst_geom
from lsst.afw.image import SKY
from lsst.afw.image import VisitInfo
from lsst.afw.image import ExposureF as ExpF
from lsst.ts.wep.bsc.CamFactory import CamFactory
from lsst.ts.wep.Utility import CamType, expandDetectorName

if __name__ == "__main__":

    # where is the telescope pointing
    ra_pointing = 0.0  # in degrees
    dec_pointing = 0.0  # in degrees
    bore_sight = lsst_geom.SpherePoint(ra_pointing*lsst_geom.degrees,
                                    dec_pointing*lsst_geom.degrees)

    rotation_angle = lsst_geom.Angle(0.0*lsst_geom.degrees)

    # must define the type of the rotation
    # taken from the code:
    #
    # SKY, ///< Position angle of focal plane +Y, measured from N through E.
    #      ///< At 0 degrees, +Y is along N and +X is along E/W depending on
    #      ///  handedness.
    #      ///< At 90 degrees, +Y is along E and +X is along S/N depending on
    #     ///   handedness.

    rotation_type = SKY

    visit_info = VisitInfo(boresightRaDec=bore_sight,
                        boresightRotAngle=rotation_angle,
                        rotType=rotation_type)

    # instantiate the LSST camera
    lsst_cam_mapper = obs_lsst.lsstCamMapper.LsstCamMapper()
    lsst_cam = lsst_cam_mapper.camera

    results = []
    results_columns = ['wcsMethod', 'conversionMethod',
                       'raft', 'sensor', 'ra', 'dec', 'xpix', 'ypix']

    for raft, sensor, detector in zip(['R22']*9,
                                      ['S00', 'S01', 'S02',
                                       'S10', 'S11', 'S12',
                                       'S20', 'S21', 'S22'],
                                      ['90', '91', '92', '93', '94',
                                       '95', '96', '97', '98']):

    # raft = 'R22'
    # sensor = 'S11'
    # detector = '94'
        print('{}_{}'.format(raft, sensor))
        #print('FlipX = {}'.format(flipX))

        wavefront_detector = lsst_cam['{}_{}'.format(raft, sensor)]

        # the flipX=True kwarg is there to get the WCS to conform to PhoSim
        # pixel coordinate conventions (I think; you may need to play around with
        # this)
        wcs = createInitialSkyWcs(visit_info, wavefront_detector, flipX=True)

        # convert a point in (RA, Dec) space to pixel coordinates
        test_ra_list = [0.0, 0.05, 0.05, -0.05, -0.05]
        test_dec_list = [0.0, 0.05, -0.05, 0.05, -0.05]
        for test_ra, test_dec in zip(test_ra_list, test_dec_list):
            ra_dec_pt = lsst_geom.SpherePoint(test_ra * lsst_geom.degrees,
                                              test_dec * lsst_geom.degrees)

            pixel_pt = wcs.skyToPixel(ra_dec_pt)
            print('Obs_LSST results Flip X')
            print('(RA, Dec)=({}, {}) to Pixel'.format(test_ra, test_dec))
            print('X: ', pixel_pt.getX())
            print('Y: ', pixel_pt.getY())
            print('')
            results.append(['obs_lsst_flipX', 'skyToPixel', raft, sensor,
                            test_ra, test_dec, pixel_pt.getX(), pixel_pt.getY()])

        # convert a point in pixel coords to RA, DEc
        xpixList = [2048, 1000, 1000, 3000, 3000]
        ypixList = [2000, 1000, 3000, 1000, 3000]
        for xpix, ypix in zip(xpixList, ypixList):
            ra_dec_pt = wcs.pixelToSky(xpix, ypix)
            print('Pixel=({}, {}) to RA, Dec'.format(xpix, ypix))
            print('RA: ', ra_dec_pt.getRa().asDegrees())
            print('Dec: ', ra_dec_pt.getDec().asDegrees())
            print('')
            results.append(['obs_lsst_flipX', 'pixelToSky', raft, sensor,
                            ra_dec_pt.getRa().asDegrees(),
                            ra_dec_pt.getDec().asDegrees(), xpix, ypix])

        wcs2 = createInitialSkyWcs(visit_info, wavefront_detector, flipX=False)

        # convert a point in (RA, Dec) space to pixel coordinates
        test_ra_list = [0.0, 0.05, 0.05, -0.05, -0.05]
        test_dec_list = [0.0, 0.05, -0.05, 0.05, -0.05]
        for test_ra, test_dec in zip(test_ra_list, test_dec_list):
            ra_dec_pt = lsst_geom.SpherePoint(test_ra * lsst_geom.degrees,
                                            test_dec * lsst_geom.degrees)

            pixel_pt = wcs2.skyToPixel(ra_dec_pt)
            print('Obs_LSST results')
            print('(RA, Dec)=({}, {}) to Pixel'.format(test_ra, test_dec))
            print('X: ', pixel_pt.getX())
            print('Y: ', pixel_pt.getY())
            print('')
            results.append(['obs_lsst', 'skyToPixel', raft, sensor,
                            test_ra, test_dec, pixel_pt.getX(), pixel_pt.getY()])

        # convert a point in pixel coords to RA, DEc
        xpixList = [2048, 1000, 1000, 3000, 3000]
        ypixList = [2000, 1000, 3000, 1000, 3000]
        for xpix, ypix in zip(xpixList, ypixList):
            ra_dec_pt = wcs2.pixelToSky(xpix, ypix)
            print('Pixel=({}, {}) to RA, Dec'.format(xpix, ypix))
            print('RA: ', ra_dec_pt.getRa().asDegrees())
            print('Dec: ', ra_dec_pt.getDec().asDegrees())
            print('')
            results.append(['obs_lsst', 'pixelToSky', raft, sensor,
                            ra_dec_pt.getRa().asDegrees(),
                            ra_dec_pt.getDec().asDegrees(), xpix, ypix])

        print('Phosim WCS Results')
        exp = ExpF('/astro/users/brycek/epyc/users/brycek/Commissioning/aos/ts_phosim/notebooks/analysis_scripts/test_output/baselineCloseLoop/input/rerun/run1/postISRCCD/09006001-g/R22/postISRCCD_09006001-g-{}-{}-det0{}.fits'.format(raft, sensor, detector))
        w = exp.getWcs()
        test_ra_list = [0.0, 0.05, 0.05, -0.05, -0.05]
        test_dec_list = [0.0, 0.05, -0.05, 0.05, -0.05]
        for test_ra, test_dec in zip(test_ra_list, test_dec_list):
            ra_dec_pt = lsst_geom.SpherePoint(test_ra * lsst_geom.degrees,
                                            test_dec * lsst_geom.degrees)

            pixel_pt = w.skyToPixel(ra_dec_pt)
            print('(RA, Dec)=({}, {}) to Pixel'.format(test_ra, test_dec))
            print('X: ', pixel_pt.getX())
            print('Y: ', pixel_pt.getY())
            print('')
            results.append(['phosim_wcs', 'skyToPixel', raft, sensor,
                            test_ra, test_dec, pixel_pt.getX(), pixel_pt.getY()])

        xpixList = [2048, 1000, 1000, 3000, 3000]
        ypixList = [2000, 1000, 3000, 1000, 3000]
        for xpix, ypix in zip(xpixList, ypixList):
            ra_dec_pt = w.pixelToSky(xpix, ypix)
            print('Pixel=({}, {}) to RA, Dec'.format(xpix, ypix))
            print('RA: ', ra_dec_pt.getRa().asDegrees())
            print('Dec: ', ra_dec_pt.getDec().asDegrees())
            print('')
            results.append(['phosim_wcs', 'pixelToSky', raft, sensor,
                            ra_dec_pt.getRa().asDegrees(),
                            ra_dec_pt.getDec().asDegrees(), xpix, ypix])

        print('sims_CoordUtils results')
        camera = CamFactory.createCam(CamType.ComCam)
        camera.setObsMetaData(ra_pointing, dec_pointing, 0.0)
        chipName = expandDetectorName('{}_{}'.format(raft, sensor))
        test_ra_list = [0.0, 0.05, 0.05, -0.05, -0.05]
        test_dec_list = [0.0, 0.05, -0.05, 0.05, -0.05]
        for test_ra, test_dec in zip(test_ra_list, test_dec_list):
            pix_x, pix_y = camera._wcs.pixelCoordsFromRaDec(test_ra, test_dec, chipName=chipName,
                                                            epoch=2000.0, includeDistortion=True)
            print('(RA, Dec)=({}, {}) to Pixel'.format(test_ra, test_dec))
            print('X: ', pix_x)
            print('Y: ', pix_y)
            print('')
            results.append(['sims_coordUtils', 'skyToPixel', raft, sensor,
                            test_ra, test_dec, pix_x, pix_y])

        xpixList = [2048, 1000, 1000, 3000, 3000]
        ypixList = [2000, 1000, 3000, 1000, 3000]
        for xpix, ypix in zip(xpixList, ypixList):
            ra_pt, dec_pt = camera._wcs.raDecFromPixelCoords(xpix, ypix, chipName,
                                                            epoch=2000.0, includeDistortion=True)
            print('Pixel=({}, {}) to RA, Dec'.format(xpix, ypix))
            print('RA: ', ra_pt)
            print('Dec: ', dec_pt)
            print('')
            results.append(['sims_coordUtils', 'pixelToSky', raft, sensor,
                            ra_pt, dec_pt, xpix, ypix])

    results_df = pd.DataFrame(results, columns=results_columns)
    results_df.to_csv('coord_convert_comcam.csv', index=False)
