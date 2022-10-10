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

import galsim
import numpy as np
from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.Utility import CamType


def getPsfGradPerZernike(
    jmax: int = 22,
    inst: Instrument = None,
) -> np.ndarray:
    """Get the gradient of the PSF FWHM with respect to each Zernike.

    Parameters
    ----------
    jmax : int
        The max Zernike Noll index, inclusive. (the default is 22.)
    inst : Instrument, optional
        The ts.wep instrument used in the Algorithm class to solve the TIE.
        (the default is the default LsstCam.)

    Returns
    -------
    np.ndarray
        Gradient of the PSF FWHM with respect to the corresponding Zernike.
        Units are arcsec / micron.
    """
    # If inst is None, get the default LsstCame
    if inst is None:
        inst = Instrument()
        inst.configFromFile(
            160,  # Donut dimension; irrelevant for this function
            CamType.LsstCam,  # Default is LsstCam
        )

    # Get the inner and outer radii of the instrument
    R_outer = inst.apertureDiameter / 2
    R_inner = R_outer * inst.obscuration

    # Calculate the conversion factors
    conversion_factors = np.zeros(jmax - 4)
    for i in range(4, jmax + 1):
        # Set coefficients for this Noll index: coefs = [0, 0, ..., 1]
        # Note the first coefficient is Noll index 0, which does not exist and
        # is therefore always ignored by galsim
        coefs = [0] * i + [1]

        # Create the Zernike polynomial with these coefficients
        Z = galsim.zernike.Zernike(coefs, R_outer=R_outer, R_inner=R_inner)

        # We can calculate the size of the PSF from the RMS of the gradient of
        # the wavefront. The gradient of the wavefront perturbs photon paths.
        # The RMS quantifies the size of the collective perturbation.
        # If we expand the wavefront gradient in another series of Zernike
        # polynomials, we can exploit the orthonormality of the Zernikes to
        # calculate the RMS from the Zernike coefficients.
        rms_tilt = np.sqrt(np.sum(Z.gradX.coef**2 + Z.gradY.coef**2))

        # Convert to arcsec per micron
        rms_tilt = np.rad2deg(rms_tilt * 1e-6) * 3600

        # Convert rms -> fwhm
        fwhm_tilt = 2 * np.sqrt(2 * np.log(2)) * rms_tilt

        # Save this conversion factor
        conversion_factors[i] = fwhm_tilt

    return conversion_factors


def convertZernikesToPsfWidth(
    zernikes: np.ndarray,
    inst: Instrument = None,
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

    Note this conversion depends on the aperture of the telescope.

    Parameters
    ----------
    zernikes: np.ndarray
        Zernike amplitudes (in microns), starting with Noll index 4.
        Either a 1D array of zernike amplitudes, or a 2D array, where each row
        corresponds to a different set of amplitudes.
    inst: Instrument, optional
        The ts.wep instrument used in the Algorithm class to solve the TIE.
        (the default is the default LsstCam.)

    Returns
    -------
    np.ndarray
        Quadrature contribution of each Zernike vector to the PSF FWHM
        (in arcseconds).
    """
    conversion_factors = getPsfGradPerZernike(zernikes.shape[-1], inst)

    # Convert the zernike amplitudes from microns to PSF FWHM
    return conversion_factors * zernikes
