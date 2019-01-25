import os
from enum import Enum
import pickle

from lsst.sims.utils import ObservationMetaData
from lsst.ts.wep.Utility import FilterType

import lsst.ts.phosim


class CamDistType(Enum):
    L1S1zer = 1
    L1S2zer = 2
    L2S1zer = 3
    L2S2zer = 4
    L3S1zer = 5
    L3S2zer = 6


class SurfaceType(Enum):
    M1 = 1
    M2 = 2
    M3 = 3
    L1F = 4
    L1B = 5
    L2F = 6
    L2B = 7
    FilterF = 8
    FilterB = 9
    L3F = 10
    L3B = 11
    FP = 12
    Chip = 13


def getWfChips():
    """Get the wavefront sensors on the focal plane.

    Returns
    -------
    list or str
        The wavefront sensor names.
    """
    return [
        ("R:0,0 S:2,2,A", "R:0,0 S:2,2,B"),
        ("R:0,4 S:2,0,A", "R:0,4 S:2,0,B"),
        ("R:4,0 S:0,2,A", "R:4,0 S:0,2,B"),
        ("R:4,4 S:0,0,A", "R:4,4 S:0,0,B")
        ]


def createObservation(obsId=None, filterType=None, boresight=None,
                       zAngleInDeg=None, rotAngInDeg=None, mjd=None):
    """Create a new Observation based on the provided parameters.

    Parameters
    ----------
    obsId : int, optional
        Observation Id. (the default is None.)
    filterType : FilterType, optional
        Active filter type. (the default is None.)
    boresight : tuple, optional
        Telescope boresight in (ra, decl). (the default is None.)
    zAngleInDeg : float, optional
        Zenith angle in degree. (the default is None.)
    rotAngInDeg : float, optional
        Camera rotation angle in degree between -90 and 90 degrees. (the
        default is None.)
    mjd : float, optional
        MJD of observation. (the default is None.)

    Returns
    _______
    ObservationMetaData
        An observation with the set parameters.
    """
    _altitude = 90
    if zAngleInDeg is not None:
        if not (0 <= zAngleInDeg <= 90):
            raise ValueError("zAngleInDeg must be between 0 and 90 degrees.")
        _altitude = 90 - zAngleInDeg
    
    _ra = 0
    _dec = 0
    if boresight is not None:
        _ra = boresight[0]
        _dec = boresight[1]

    _filterString = FilterType.REF.toString()
    if filterType is not None:
        _filterString = filterType.toString()

    _rotAngInDeg = 0
    if rotAngInDeg is not None:
        _rotAngInDeg = rotAngInDeg

    _obsId = 9999
    if obsId is not None:
        _obsId = obsId

    _mjd = 59552.3
    if mjd is not None:
        _mjd = mjd

    obs = ObservationMetaData(
            mjd=_mjd, 
            pointingRA=_ra,
            pointingDec=_dec,
            rotSkyPos=_rotAngInDeg,
            bandpassName=_filterString
        )
    obs.OpsimMetaData = {
        "obsHistID": _obsId,
        "altitude": _altitude
    }
    return obs


def getOpsimObservation(tract=0, target=0):
    """Get the Opsim Observation corresponding to the tract and target.
    The tracts are shown in doc/tracts.png.

    Parameters
    ----------
    tract : int, optional
        The tract to get the catalog from. 
        (the default is 0. See doc/tracts.png)
    target : int, optional
        The target within a tract to get the catalog from. 
        (the default is 0. See doc/tracts.png)

    Returns
    -------
    ObservationMetaData
        The selected observation.
    """
    tractFile = os.path.join(getModulePath(), "configData/tract/{}".format(tract))
    with open(tractFile, "rb") as r:
        tract = pickle.load(r, encoding="latin1")
    return tract[target]


def getModulePath(module=lsst.ts.phosim, startIdx=1, endIdx=-4):
    """Get the path of module.

    Parameters
    ----------
    module : str, optional
        Module name. (the default is lsst.ts.ofc.)
    startIdx : int, optional
        Start index. (the default is 1.)
    endIdx : int, optional
        End index. (the default is -4.)

    Returns
    -------
    str
        Directory path of module based on the start and end indexes.
    """

    # Get the path of module
    modulePathList = os.path.dirname(module.__file__).split(
                                os.sep)[int(startIdx):int(endIdx)]
    modulePath = os.path.join(os.sep, *modulePathList)

    return modulePath


def opt2ZemaxCoorTrans(x, y, z):
    """Do the coordinate transformation from optics to Zemax.

    Parameters
    ----------
    x : float
        x coordinate in optics.
    y : float
        y coordinate in optics.
    z : float
        z coordinate in optics.

    Returns
    -------
    float
        x coordinate in Zemax.
    float
        y coordinate in Zemax.
    float
        z coordinate in Zemax.
    """

    return -x, y, -z


def zemax2optCoorTrans(x, y, z):
    """Do the coordinate transformation from Zemax to optics.

    Parameters
    ----------
    x : float
        x coordinate in Zemax.
    y : float
        y coordinate in Zemax.
    z : float
        z coordinate in Zemax.

    Returns
    -------
    float
        x coordinate in optics.
    float
        y coordinate in optics.
    float
        z coordinate in optics.
    """

    return -x, y, -z


def mapSurfNameToEnum(surfName):
    """Map the surface name to Enum: SurfaceType. 

    Parameters
    ----------
    surfName : str
        Surface name.

    Returns
    -------
    SurfaceType
        SurfaceType enum.

    Raises
    ------
    ValueError
        No surface type is found.
    """

    nameList = [surfaceType.name for surfaceType in SurfaceType]
    try:
        idx = nameList.index(surfName.strip())
        return SurfaceType(idx + 1)
    except ValueError:
        raise ValueError ("No surface type is found.")


if __name__ == "__main__":
    pass
