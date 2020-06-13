# common functions for AOS analysis
import os
import lsst.daf.persistence as dafPersist
from astropy.table import Table
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.cm as cm 

def readPostISRImage(data_dir, focalType = 'extra', obsId=None, raft = 'R22',
                     detector = 'S00', detNum = None, verbose=True):
    ''' A function to read the post ISR image for a given CCD (sensor) 
    using Butler (so it has to be a butler repository). 

    Parameters:
    -----------
    data_dir : the location of postISR data: a path that contains 
               /input/rerun/run1/postISRCCD/
    focalType: 'extra' or 'intra' 
    obsId: id corresponding to the focal type, if comcam, then read from a dictionary:
         {'intra':9006001,  'extra':9006002} , otherwise set manually 
    detector: ccd name, one of ['S00', 'S01', 'S02','S10', 'S11', 'S12', 'S20', 'S21', 'S22']. 
         'S00' by default 
    raft : str, by default 'R22'  (comcam)
    detNum : int, detector Id - by default it's read from a following dictionary for comCam
             {'S00':90, 'S01':91, 'S02':92, 'S10':93, 'S11':94, 'S12':95, 'S20':96, 'S21':97, 'S22':98}
    verbose : boolean, True by default  - whether to print info about the progress of reading images 
    
    Returns:
    --------
    image :  an NxM array with the CCD postISR image 

    '''
    # Read in the postISR image using the Butler 
    obsIdDic = {'intra':9006001,  'extra':9006002} 
    
    if not obsId: # if not provided, reading it from a dict
        obsId = obsIdDic[focalType]
    
    detectors = ['S00', 'S01', 'S02','S10', 'S11', 'S12', 'S20', 'S21', 'S22']
    sensor = raft+'_'+detector 

    detNumDict = {'S00':90, 'S01':91, 'S02':92, 'S10':93, 'S11':94, 
                  'S12':95, 'S20':96, 'S21':97, 'S22':98}
    if not detNum:
        detNum = detNumDict[detector]

    data_id = {'visit': obsId, 'filter': 'g', 'raftName': raft, 
               'detectorName': detector, 'detector': detNum
              }


    # Read each figure as a postage stamp, store data to an array 
    if verbose: 
        print('\nReading data from %s'%data_dir)
        print('For sensor %s '%sensor)
    repo_dir = os.path.join(data_dir, 'input/rerun/run1')
    butler = dafPersist.Butler(repo_dir)

    # show what keys are needed by the `postISRCCD` data type.... 
    # butler.getKeys('postISRCCD')
    # yields {'visit': int, 'filter': str,'raftName': str, 'detectorName': str, 'detector': int}
    post = butler.get('postISRCCD', **data_id) 

    # store in a dictionary
    image = post.image.array

    if verbose: print('Done\n')
    
    return image 






def readCentroidInfo(data_dir, focalType='extra', raft='R22',detector='S00'):
    ''' Read the centroid file info 

    Parameters:
    ------------
    data_dir : the location of postISR data: a path that contains 
                   /input/rerun/run1/postISRCCD/
    focalType: 'extra' or 'intra' 
    
    Returns:
    ---------
    centroids : an astroPy Table with a list of centroids 
    centFlag:  a flag that is False if the file is not found , True  otherwise
    '''
    # read centroid info 
    centr_dir = os.path.join(data_dir, 'iter0','img',focalType)
    print('Reading centroid files from %s'%centr_dir)
    fname = 'centroid_lsst_e_%d_f1_%s_%s_E000.txt'%(obsId,raft,detector)
    path = os.path.join(centr_dir,fname)
    try:
        print('Using  %s'%fname)
        centroid = Table.read(centr_dir+'/'+fname, format='ascii')
        centFlag = True
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        centFlag = False
        centroid = []
    return centroid, centFlag    

def readPostageStars(data_dir, ):
    # read postage stars info
    # Read in the postage image catalog
    try:
        fname = 'postagedonutStarsExtraIntra.txt'
        postage = Table.read(os.path.join(data_dir,fname), format='ascii')
        print('Reading info about postage-stamp images from %s'%fname)
        postFlag = True
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        postFlag = False
    return postage, postFlag

# helper funtion for the colorbar 
# from https://joseph-long.com/writing/colorbars/
def colorbar(mappable):
    last_axes = plt.gca()
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar


def plotImage(image,ax=None, log=False, sensor='R22_S00', focalType='extra',
             postage=None,postFlag=False, centroid=None, centFlag=False, 
              Nstars=2,starMarker='redCross',starMarkerArgs=None,
             centMarkerArgs = None,centMarker='redCross',
             starLabelArgs=None):
    ''' A function  to plot a CCD image

    Parameters:
    -----------
    image : an NxM image array with pixel values . 
        It will be transposed for plotting 
    ax : plt.axis,  an axis to plot the image on
    log : bool, False by default - whether to pot log(counts) or not 
    

    Returns:
    --------
    None

    '''
    
    if not ax:
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        
    # plot the image 
    if log : 
        
        plottable = np.log10(image.T)
        cbar_label = r'$\log_{10}(\mathrm{counts})$'
    else:
        plottable = image.T
        cbar_label = r'$\mathrm{counts}$'
    
        img = ax.imshow(plottable,# vmin = 2.45, vmax=2.75,
                  origin='lower')
        cbar= colorbar(mappable=img)
        cbar.set_label(label=cbar_label, weight='normal', )
        ax.set_xlabel('x [px]')
        ax.set_ylabel('y [px]')
        ax.set_title('postISR image, sensor %s, %s-focal'%(sensor,focalType))

    # try to figure out how many stars there are from postage file 
    if postFlag: 
        m1 = postage['abbrevDetectorName'] == sensor
        m2 = postage['focalPlane'] == focalType
        mask = m1 * m2 
        Nstars = len(postage[mask])

    starMarkerArgDict = {'yellowSquare': {'marker':'s', 'markersize':40, 'markerfacecolor':'none',
                                          'markeredgewidth':2.5, 'markeredgecolor':'y',
                                         } ,
                        'redCross':{'marker':'+', 'markersize':10,'markeredgecolor':'red', 
                                    'markerfacecolor':'red', 'markeredgewidth':2,
                                   }
                        }
    if starMarkerArgs is None:
        starMarkerArgs = starMarkerArgDict[starMarker]
    if centMarkerArgs is None:
        centMarkerArgs = starMarkerArgDict[centMarker]
    if starLabelArgs  is None:
        starLabelArgs = {'fontsize':16, 'color':'white'}
        
    for i in range(Nstars):
        if postFlag:
            x,y = postage[mask]['xpos'][i], postage[mask]['ypos'][i]
            starId = postage[mask]['starId'][i]
            #print(starId,x,y)
            ax.plot(x,y,  **starMarkerArgs)
            ax.text(x-40, y-180, starId, **starLabelArgs)

        if centFlag:
            #get the centroid location for that star
            starCentroid = centroid[centroid['SourceID'].data.astype(int) == starId]
            xCen, yCen = starCentroid['AvgX'].data  , starCentroid['AvgY'].data 
            
            #plot it with a cross 
            ax.plot(xCen, yCen, **centMarkerArgs)

    plt.tight_layout()

    
    
    
def plotZernikesAndCCD(image, rmsErrors, sepInPerc=10, xlims=[1525,2025], ylims=[750,1250],
                      sensor = 'R22_S00'):
    '''  Function to plot both rms zernike errors and CCD image as a two-panel plot,
    restricting the CCD image to show only the relevant donut 
    
    Parameters:
    -----------
    image: NxM array with CCD data to plot 
    rmsErrors : an array with size 
    
    Returns:
    -------
    None
    '''
    xmin,xmax = xlims[0], xlims[1]
    ymin,ymax = ylims[0], ylims[1]

    fig, ax = plt.subplots(1,2,figsize=(16,6))
    
    sepInRadii = sepInPercToRadii(sepInPerc)
    
    
    if np.shape(rmsErrors)[0] == 19 : 
        ax[0].plot(np.arange(19)+4, rmsErrors, 
         '-o', lw=3, )# color = cmap(colors[i]))
    else:
        print('Need the Zernike rms errors to be an array with 19 elements')

    ax[0].set_xlabel('Zernike Number', size=18)
    ax[0].set_ylabel('RMS WFS vs OPD (microns)', size=18)
    ax[0].set_title('%s, Star Sep = %.1f radii' % (sensor,sepInRadii), 
                    size=18)

    # plot the postage stamp
    img = ax[1].imshow(np.log10(image[ymin:ymax, xmin:xmax]), vmin = 0.01,
                       cmap=cm.get_cmap('Greys'),
              origin='lower')
    #cbar= colorbar(mappable=img)
    #cbar.set_label(label=r'$\log_{10}(\mathrm{counts})$', weight='normal', )
    ax[1].set_xlabel('x [px]')
    ax[1].set_ylabel('y [px]')
    ax[1].set_title('postISR image')
    plt.tight_layout()
    plt.savefig('img_AOS_singleAmpSep_postIsr_sep%d.png'%sepInPerc,
                bbox_inches='tight', dpi=150)

    
def sepInPercToRadii(sepInPerc):
    ''' Function to convert separation in percentage of amplifier ra span
    to donut radius 
    
    Parameters:
    ----------
    sepInPerc : int - separation in percent (eg. 10)
    
    Returns:
    ---------
    sepInRadii : float  - separation in donut radii
    '''
                    
    yPxAmp = 2048 # in pixels
    donutPxRadius = 66 # px            
    sepInPx = sepInPerc*0.01*yPxAmp 
    sepInRadii = sepInPx / donutPxRadius
                    
    return sepInRadii             


def plotPostageStamps(data_dir, sensor='R22_S00', focalType='extra', Nstars=2,
                     cbarX0Y0DxDy = [0.13, 0.06, 0.76, 0.01],sepInPerc=3
                     ):

    sepInRadii = sepInPercToRadii(sepInPerc)
    print(sepInRadii)
    
    imgType = ['singleSciImg','imgDeblend_full', 'imgDeblend_resized']
    postage_dir = os.path.join(data_dir, 'postage')
    print('Using postage images from %s'%postage_dir)
    suptitle = '%s %s-focal, sep=%.1f donut radii'%(sensor,focalType,sepInRadii)
     
    fig,ax = plt.subplots(Nstars,len(imgType),figsize=(12,4*Nstars))

    for col in range(len(imgType)): # columns : each imgType is one column 
        ax[0,col].set_title(imgType[col], fontsize=17)
        for row in range(Nstars): # Nstars   rows : one per star
            fname = focalType+'_'+imgType[col]+"_sensor-"+sensor+"_star-"+str(row)+'_'
            #print(fname)
            for x in os.listdir(postage_dir):
                if x.startswith(fname): 
                    print(x)
                    fname = x
                    word = 'id'
                    loc = fname.find(word)
                    starId = fname[loc+len(word)+1]
            image = np.loadtxt(postage_dir+'/'+fname)
            if image.ndim == 2  :
                mappable = ax[row,col].imshow(image)
                ax[row,col].text(0.1,0.1,'star %d, id %s'%(row,starId) , 
                                 fontsize=17, color='white', 
                                 transform=ax[row,col].transAxes)
            else:
                ax[row,col].remove()
                ax[row,col].text(0.2,0.5, 'image.ndim < 2 ',fontsize=15,
                            transform=ax[row,col].transAxes)

    # that's for horizontal cbar on the bottom 
    cbar_ax = fig.add_axes(cbarX0Y0DxDy)     #  (x0 ,y0  , dx,  dy )  
    cbar = fig.colorbar(mappable, cax = cbar_ax,  orientation='horizontal')                    
    cbar.set_label(label='counts',weight='normal', fontsize=17)

    fig.suptitle(suptitle, fontsize=17)
    plt.savefig('img_AOS_singleAmpSep_'+sensor+'_'+focalType+'_postageImg.png', 
                bbox_inches='tight', dpi=100)