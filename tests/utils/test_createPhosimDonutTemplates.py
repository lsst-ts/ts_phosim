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

import unittest
import os
import shutil
import tempfile
import numpy as np

from lsst.ts.phosim.utils.Utility import getConfigDir, getModulePath
from lsst.ts.phosim.utils.CreatePhosimDonutTemplates import CreatePhosimDonutTemplates
from lsst.ts.wep.utility import DefocalType


class TestCreatePhosimDonutTemplates(unittest.TestCase):
    """Test the CreatePhosimDonutTemplates class."""

    def setUp(self):
        modulePath = getModulePath()
        self.testDataDir = os.path.join(modulePath, "tests", "testData")

        # Location where phosim input files exist.
        # Also where temp work directory is created.
        self.templateDataDir = os.path.join(getConfigDir(), "donutTemplateData")

        # Location where templates are created
        self.templateDir = os.path.join(self.templateDataDir, "phosimTemplates")

        # Temporary work directory inside ts_wep/tests/
        testDir = os.path.join(modulePath, "tests")
        self.tempTestDirectory = tempfile.TemporaryDirectory(dir=testDir)
        self.tempWorkDir = self.tempTestDirectory.name
        self.createPhosimDonuts = CreatePhosimDonutTemplates(self.tempWorkDir)

    def tearDown(self):
        # Remove test template from template directory
        try:
            os.remove(os.path.join(self.templateDir, "extra_template-R22_S10.txt"))
        except FileNotFoundError:
            pass

        # Remove work directory in case of clean up method test failure
        self.tempTestDirectory.cleanup()

    def _copyPhosimFiles(self):
        # copy the raw amp image file
        shutil.copy(
            os.path.join(
                self.testDataDir,
                "testDonutTemplates",
                "phosimOutput",
                "MC_H_20211231_006001_R22_S10.fits",
            ),
            os.path.join(self.tempWorkDir, "raw", "MC_H_20211231_006001_R22_S10.fits"),
        )

        # copy the centroid file ( for a different visit )
        shutil.copy(
            os.path.join(
                self.testDataDir,
                "testDonutTemplates",
                "phosimOutput",
                "centroid_lsst_e_9006001_f1_R22_S10_E000.txt",
            ),
            os.path.join(self.tempWorkDir, "phosimOutput", "extra"),
        )

    def testCreateAndCleanUpWorkDirectories(self):
        # Test clean up of work directories
        self.createPhosimDonuts.cleanUpWorkDirs()
        self.assertFalse(os.path.exists(self.tempWorkDir))

        # Test creation of work directories
        self.createPhosimDonuts.createWorkDirectories()
        self.assertTrue(os.path.exists(self.tempWorkDir))
        self.assertTrue(
            os.path.exists(os.path.join(self.tempWorkDir, "phosimOutput", "extra"))
        )

    def testCreateDetectorLists(self):
        testDetectors = "R22_S00 R22_S11"

        detStrPhosim = self.createPhosimDonuts.createDetectorLists(testDetectors)

        self.assertEqual(detStrPhosim, "R22_S00|R22_S11")

    def testIngestImages(self):
        # Populate the raw phosim output directory
        self.createPhosimDonuts.createWorkDirectories()
        self._copyPhosimFiles()
        # Run the ingest
        self.createPhosimDonuts.ingestImages()
        registry = self.createPhosimDonuts.butler.registry

        print(list(registry.queryDataIds('exposure')))
        print(list(registry.queryCollections()))
        self.assertEqual(
            list(registry.queryDataIds("exposure"))[0]["exposure"], 4021123106001
        )

    def testRunISR(self):
        self.createPhosimDonuts.createWorkDirectories()
        # Populate the repo
        self._copyPhosimFiles()
        self.createPhosimDonuts.ingestImages()

        # Run the ISR
        self.createPhosimDonuts.runISR("")

        # Query registry
        registry = self.createPhosimDonuts.butler.registry
        datasetsList = list(
            registry.queryDatasets(
                "postISRCCD",
                collections=[f"{self.createPhosimDonuts.outputRunName}"],
            )
        )

        self.assertEqual(len(datasetsList), 1)

    def testGenerateDataIdLists(self):
        self.createPhosimDonuts.createWorkDirectories()
        # Populate the repo
        self._copyPhosimFiles()
        self.createPhosimDonuts.ingestImages()
        # Run the ISR
        self.createPhosimDonuts.runISR("")

        intraIdList, extraIdList = self.createPhosimDonuts.generateDataIdLists(
            9006002, 9006001
        )

        self.assertEqual(intraIdList, [])
        self.assertCountEqual(
            extraIdList,
            [
                {
                    "band": "g",
                    "instrument": "LSSTCam",
                    "detector": 93,
                    "physical_filter": "g_6",
                    "exposure": 4021123106001,
                }
            ],
        )

    def testCutOutTemplatesAndSave(self):
        # Move centroid file into place
        self.createPhosimDonuts.createWorkDirectories()
        self._copyPhosimFiles()
        self.createPhosimDonuts.ingestImages()
        # Run the ISR
        self.createPhosimDonuts.runISR("")

        intraIdList, extraIdList = self.createPhosimDonuts.generateDataIdLists(
            9006002, 9006001
        )
        self.createPhosimDonuts.cutOutTemplatesAndSave(
            240,
            DefocalType.Extra,
            extraIdList,
        )

        newTemplate = np.genfromtxt(
            os.path.join(self.templateDir, "extra_template-R22_S10.txt")
        )
        trueTemplate = np.genfromtxt(
            os.path.join(
                self.testDataDir, "testDonutTemplates", "extra_template-R22_S10.txt"
            )
        )
        np.testing.assert_array_equal(newTemplate, trueTemplate)


if __name__ == "__main__":
    # Do the unit test
    unittest.main()
