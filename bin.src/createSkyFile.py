#!/usr/bin/env python

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

import argparse
import numpy as np
import pandas as pd
from astropy import units as u

from lsst.daf.butler import Butler

from lsst.ts.wep.task.RefCatalogInterface import RefCatalogInterface


class CreateSkyFile(object):
    """
    Code to create input sky catalog files for the closed loop
    from reference catalogs available in a butler repository.

    Parameters
    ----------
    butlerPath : str
        Filepath to the butler repository.
    collectionName : str
        Collection name inside the repository that contains
        the reference catalogs.
    refCatName : str
        The name of the reference catalog dataset.
    """

    def __init__(self, butlerPath, collectionName, refCatName):
        self.butler = Butler(butlerPath)
        self.collection = collectionName
        self.refCatName = refCatName

    def createSkyFile(
        self, ra, dec, filterName, outFile, radius=1.8, magMax=99.0, magMin=-99.0
    ):
        """
        Get a catalog of reference catalog sources that cover
        a given part of the sky and save it in a file that can
        be used as input for the closed loop.

        Parameters
        ----------
        ra : float
            RA of pointing center in degrees.
        dec : float
            Dec of pointing center in degrees.
        filterName : str
            Name of the filter to simulate in the closed loop.
        outFile : str
            Name of the output catalog file.
        radius : float, optional
            Radius in degrees of the sky footprint to write out to file.
            (the default is 1.8).
        magMax : float, optional
            Maximum magnitude of sources to include in the closed
            loop catalog. (the default is 99.0).
        magMin : float, optional
            Minimum magnitude of sources to include in the closed
            loop catalog. (the default is -99.0).
        """
        refCatInterface = RefCatalogInterface(ra, dec, 0.0)
        shardIds = refCatInterface.getShardIds(radius=radius)
        dataRefs, dataIds = refCatInterface.getDataRefs(
            shardIds, self.butler, self.refCatName, self.collection
        )
        catColumns = ["Ra", "Decl", "Mag"]
        skyCat = pd.DataFrame([], columns=catColumns)
        for dataRef in dataRefs:
            shardCat = self.butler.getDirect(dataRef.ref).asAstropy().to_pandas()
            filterFlux = shardCat[f"{filterName}_flux"].values * u.nJy
            filterMag = filterFlux.to(u.ABmag)
            raDeg = np.rad2deg(shardCat["coord_ra"])
            decDeg = np.rad2deg(shardCat["coord_dec"])
            subCat = pd.DataFrame(
                (np.array([raDeg, decDeg, filterMag]).T), columns=catColumns
            )
            if skyCat is None:
                skyCat = subCat
            else:
                skyCat = pd.concat([skyCat, subCat], ignore_index=True)
        else:
            skyCat = skyCat.query(f"Mag < {magMax} and Mag > {magMin}")

        skyCat.to_csv(outFile, index_label="# Id", sep="\t")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("butlerPath", type=str, help="Butler Repository with refcats.")
    parser.add_argument(
        "collectionName", type=str, help="Collection with refcats within butler repo."
    )
    parser.add_argument("refCatName", type=str, help="Reference Catalog Dataset Name.")
    parser.add_argument("ra", type=float, help="RA of pointing in degrees.")
    parser.add_argument("dec", type=float, help="Dec of pointing in degrees.")
    parser.add_argument(
        "filter", type=str, help="Filter name for magnitudes from source catalog."
    )
    parser.add_argument("outFile", type=str, help="Output filepath for final catalog.")
    parser.add_argument(
        "--radius", type=float, default=1.8, help="Radius of catalog footprint."
    )
    parser.add_argument(
        "--magMax",
        type=float,
        default=99.0,
        help="Maximum magnitude of sources in catalog.",
    )
    parser.add_argument(
        "--magMin",
        type=float,
        default=-99.0,
        help="Minimum magnitude of sources in catalog.",
    )

    args = parser.parse_args()
    skyCreator = CreateSkyFile(args.butlerPath, args.collectionName, args.refCatName)
    skyCreator.createSkyFile(
        args.ra,
        args.dec,
        args.filter,
        args.outFile,
        radius=args.radius,
        magMax=args.magMax,
        magMin=args.magMin,
    )
