import os
from enum import Enum

from lsst.utils import getPackageDir


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


def getModulePath():
    """Get the path of module.

    Returns
    -------
    str
        Directory path of module.
    """

    return getPackageDir("ts_phosim")


def getConfigDir():
    """Get the directory of configuration files.

    Returns
    -------
    str
        Directory of configuration files.
    """

    return os.path.join(getModulePath(), "policy")


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
    enum 'SurfaceType'
        SurfaceType.

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
        raise ValueError("No surface type is found.")


if __name__ == "__main__":
    pass
