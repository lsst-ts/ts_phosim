# This file is part of ts_phosim.
#
# Developed for the LSST Telescope and Site Systems.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os

import numpy as np
from lsst.ts.phosim.utils.Utility import getConfigDir
from lsst.ts.wep.ParamReader import ParamReader
from lsst.ts.wep.Utility import CamType


def _load_conversion_factors() -> dict:
    """Load the conversion factors for ConvertZernikesToPsfWidth.

    Conversion factors saved in
    `ts_phosim/policy/convertZernikesToPsfWidthSetting.yaml`
    """
    # load the settings file
    settingFilePath = os.path.join(
        getConfigDir(), "convertZernikesToPsfWidthSetting.yaml"
    )
    settingFile = ParamReader(filePath=settingFilePath)

    # get all conversion factors
    conversion_factors = settingFile.getSetting("conversion_factors")
    conversion_factors = {
        getattr(CamType, camType): np.array(factors)
        for camType, factors in conversion_factors.items()
    }

    return conversion_factors


def convertZernikesToPsfWidth(
    zernikes: np.ndarray,
    camType: CamType = CamType.LsstCam,
) -> np.ndarray:
    """Convert Zernike amplitudes to quadrature contribution to the PSF FWHM.

    Converting Zernike amplitudes to their quadrature contributions to the PSF
    FWHM allows for easier physical interpretation of Zernike amplitudes and
    the performance of the AOS system.

    For example, if ConvertZernikesToPsfWidth([Z4, Z5, Z6]) = [0.1, -0.2, 0.3],
    then these Zernike perturbations increase the PSF FWHM by
    sqrt[(0.1)^2 + (-0.2)^2 + (0.3)^2] ~ 0.37 arcsecs.

    If the AOS perfectly corrects for these perturbations, the PSF FWHM will
    not increase in size. However, imagine that the AOS estimates and corrects
    for ConvertZernikesToPsfWidth([Z4, Z5, Z6]) = [0.1, -0.3, 0.2].
    Then the resultant PSF FWHM will be degraded by
    sqrt[(0.1 - 0.1)^2 + (-0.2 - (-0.3))^2 + (0.3 - 0.2)^2] ~ 0.14 arcsecs.

    Note that this conversion is camera dependent. Currently, conversions exist
    for LsstCam (the main Rubin camera) and AuxTel (the camera on the Auxiliary
    Telescope).

    Parameters
    ----------
    zernikes: np.ndarray
        Zernike amplitudes (in microns), starting with Noll index 4. Can
        include up to Noll index 37 (inclusive).
    camType: enum 'CamType', default=CamType.LsstCam
        Camera for which the conversion is performed. Currently the valid
        options are CamType.LsstCam and CamType.AuxTel.
        Specified via lsst.ts.wep.Utility.CamType.

    Returns
    -------
    np.ndarray
        Quadrature contribution of each Zernike vector to the PSF FWHM
        (in arcseconds).

    Raises
    ------
    KeyError
        If the camType is not one of the supported cameras.
    ValueError
        If the length of the zernike array is greater than 34, since only Noll
        indices 4 -> 37 (inclusive) are supported.
    """
    # load all the conversion factors
    conversion_factors = _load_conversion_factors()

    # pull out the conversion factors for the requested camera
    try:
        conversion_factors = conversion_factors[camType]
    except KeyError as error:
        raise KeyError(
            "Currently, conversions are only available for "
            f"{', '.join([cam.name for cam in conversion_factors])}."
        ) from error

    # make sure zernike array isn't too long
    if len(zernikes) > len(conversion_factors):
        raise ValueError(
            f"Zernike conversion for {camType.name} is only supported for "
            f"Noll indices 4 -> {4 + len(conversion_factors) - 1}."
        )

    # convert the zernike amplitudes from microns
    # to PSF FWHM quadrature contributions in arcsecs
    psf_fwhm_contrib = zernikes * conversion_factors[: len(zernikes)]

    return psf_fwhm_contrib
