import numpy as np
import matplotlib
# Must be before importing matplotlib.pyplot or pylab!
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plotResMap(zfInMm, xfInMm, yfInMm, outerRinMm, resFile=None,
               writeToResMapFilePath=None):
    """Plot the mirror residue map.

    Parameters
    ----------
    zfInMm : numpy.ndarray
        Surface map in mm.
    xfInMm : numpy.ndarray
        X position in mm.
    yfInMm : numpy.ndarray
        Y position in mm.
    outerRinMm : float
        Outer radius in mm.
    resFile : str, optional
        File path of the grid surface residue map. (the default is None.)
    writeToResMapFilePath : str, optional
        File path to save the residue map. (the default is None.)
    """

    # Plot the figure
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))

    # The input data to gridSamp.m is in mm (zemax default)
    sc = ax[1].scatter(xfInMm, yfInMm, s=25, c=zfInMm*1e6, marker=".",
                       edgecolor="none")
    ax[1].axis("equal")
    ax[1].set_title("Surface map on FEA grid (nm)")
    ax[1].set_xlim([-outerRinMm, outerRinMm])
    ax[1].set_ylim([-outerRinMm, outerRinMm])
    ax[1].set_xlabel("x (mm)")

    if (resFile is not None):

        # Get the data
        data = np.loadtxt(resFile)
        NUM_X_PIXELS, NUM_Y_PIXELS, delxInMm, delyInMm = data[0, :]
        NUM_X_PIXELS = int(NUM_X_PIXELS)
        NUM_Y_PIXELS = int(NUM_Y_PIXELS)

        # Get the zp data
        zpTemp = data[1:, 0]
        zp = np.zeros((NUM_X_PIXELS, NUM_Y_PIXELS))
        for jj in range(1, NUM_X_PIXELS + 1):
            for ii in range(1, NUM_Y_PIXELS + 1):
                zp[NUM_X_PIXELS+1-jj-1, ii-1] = zpTemp[(jj-1)*NUM_X_PIXELS +
                                                       (ii-1)]

        # Minimum x and y
        minx = -0.5*(NUM_X_PIXELS-1)*delxInMm
        miny = -0.5*(NUM_Y_PIXELS-1)*delyInMm

        xx = np.linspace(minx, -minx, NUM_X_PIXELS)
        yy = np.linspace(miny, -miny, NUM_Y_PIXELS)
        xp, yp = np.meshgrid(xx, yy)

        xp = xp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS, 1))
        yp = yp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS, 1))
        zp = zp.reshape((NUM_X_PIXELS*NUM_Y_PIXELS, 1))

        sc = ax[0].scatter(xp, yp, s=25, c=zp*1e6, marker=".",
                           edgecolor="none")

    ax[0].axis("equal")
    ax[0].set_title("grid input to ZEMAX (nm)")
    ax[0].set_xlim([-outerRinMm, outerRinMm])
    ax[0].set_ylim([-outerRinMm, outerRinMm])
    ax[0].set_xlabel("x (mm)")
    ax[0].set_ylabel("y (mm)")

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    fig.colorbar(sc, cax=cbar_ax)

    if (writeToResMapFilePath is not None):
        plt.savefig(writeToResMapFilePath)
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    pass
