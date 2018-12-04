import os
from enum import Enum

import lsst.ts.phosim


class CamDistType(Enum):
    L1S1zer = 1
    L1S2zer = 2
    L2S1zer = 3
    L2S2zer = 4
    L3S1zer = 5
    L3S2zer = 6


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


def phosim2ZemaxCoorTrans(x, y, z):
    """Do the coordinate transformation from PhoSim to Zemax.

    Parameters
    ----------
    x : float
        x coordinate in PhoSim.
    y : float
        y coordinate in PhoSim.
    z : float
        z coordinate in PhoSim.

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


def zemax2phosimCoorTrans(x, y, z):
    """Do the coordinate transformation from Zemax to PhoSim.

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
        x coordinate in PhoSim.
    float
        y coordinate in PhoSim.
    float
        z coordinate in PhoSim.
    """

    return -x, y, -z


if __name__ == "__main__":
    pass
