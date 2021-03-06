.. py:currentmodule:: lsst.ts.phosim

.. _lsst.ts.phosim-version_history:

##################
Version History
##################

.. _lsst.ts.phosim-1.2.9:

-------------
1.2.9
-------------

* Use the latest **ts_wep** that removes the dependency of ``sims`` package.
* Add the Deprecation warning to unused arguments of ``epoch``, ``includeDistortion``, and ``mjd`` in **SkySim**: ``addStarByChipPos()`` and ``setObservationMetaData()``.
* Fix the scripts in ``examples/`` directory.

.. _lsst.ts.phosim-1.2.8:

-------------
1.2.8
-------------

* Remove the dependency of ``sims`` package by letting the **SkySim** class to depend on **WcsSol** class in **ts_wep**.

.. _lsst.ts.phosim-1.2.7:

-------------
1.2.7
-------------

* Use the ``sims_w_2020_38``.
* Replace the **comcamCloseLoop.py** with the **imgCloseLoop.py**.
* Update the class diagram.
* Deprecation warning:

1. Use ``setWgtAndFieldXyOfGQ()`` to replace ``setDefaultLsstGQ()`` and ``setDefaultComcamGQ()`` in **OpdMetrology.py**.
2. Use ``getOpdArgsAndFilesForPhoSim()`` to replace ``getComCamOpdArgsAndFilesForPhoSim()`` in **PhosimCmpt.py**.
3. Use ``getPistonCamStarArgsAndFilesForPhoSim()`` to replace ``getComCamStarArgsAndFilesForPhoSim()`` in **PhosimCmpt.py**.
4. Use ``analyzeOpdData()`` to replace ``analyzeComCamOpdData()`` in **PhosimCmpt.py**.
5. Use ``repackagePistonCamImgs()`` to replace ``repackageComCamAmpImgFromPhoSim()`` and ``repackageComCamEimgFromPhoSim()`` in **PhosimCmpt.py**.

.. _lsst.ts.phosim-1.2.6:

-------------
1.2.6
-------------

* Add the **CloseLoopTask** class.

.. _lsst.ts.phosim-1.2.5:

-------------
1.2.5
-------------

* Use the ``sims_w_2020_36``.

.. _lsst.ts.phosim-1.2.4:

-------------
1.2.4
-------------

* Use the ``sims_w_2020_28``.
* Removed the unused force files.

.. _lsst.ts.phosim-1.2.3:

-------------
1.2.3
-------------

* Reformat the code by ``black``.
* Add the ``black`` check to ``.githooks``.
* Ignore ``flake8`` check of E203 ans W503 for the ``black``.
* Use the ``sims_w_2020_21``.

.. _lsst.ts.phosim-1.2.2:

-------------
1.2.2
-------------

* Use ``sims_w_2020_15``.
* Use the update bending mode and grid files of M1M3 and M2.
* Update the M2 FEA correction (gravity and temperature) for the fitting of x, y coordinate in grid file.

.. _lsst.ts.phosim-1.2.1:

-------------
1.2.1
-------------

* Use ``sims_w_2020_14``.

.. _lsst.ts.phosim-1.2.0:

-------------
1.2.0
-------------

* Use ``sims_w_2020_04``.

.. _lsst.ts.phosim-1.1.9:

-------------
1.1.9
-------------

* Use ``sims_w_2019_50``.

.. _lsst.ts.phosim-1.1.8:

-------------
1.1.8
-------------

* Use ``sims_w_2019_38``.

.. _lsst.ts.phosim-1.1.7:

-------------
1.1.7
-------------

* Use ``sims_w_2019_31``.
* Use the latest versions of **ts_wep** and **ts_ofc**.
* Remove the ``conda`` package installation in **Jenkinsfile**.
* Update the permission of workspace after the unit test.

.. _lsst.ts.phosim-1.1.6:

-------------
1.1.6
-------------

* Use ``sims_w_2019_29``.
* Supress the warning in unit tests.
* Fix the warning of nan in atmosphere structure function.
* Rotate the OPD and support the sky file, minimum DOF, and M1M3 force error ratio in command line tasks.

.. _lsst.ts.phosim-1.1.5:

-------------
1.1.5
-------------

* Use ``sims_w_2019_24``.
* Support the eimage in **comcamCloseLoop.py**.
* Depend on the **SensorWavefrontError** in **ts_wep**.
* Update the table file.

.. _lsst.ts.phosim-1.1.4:

-------------
1.1.4
-------------

* Minor bugs fixed.
* Add the get methods for **SkySim** and **OpdMetrology** classes.
* Use the **CamType** of **ts_wep** module in **TeleFacade** class.
* Update **PhosimCmpt** class to use the interface classes of **ts_wep** and **ts_ofc**.
* Use the scientific pipeline of ``sims_w_2019_20``.
* Add the command line tasks of close-loop simulation.

.. _lsst.ts.phosim-1.1.3:

-------------
1.1.3
-------------

* Combine with **ts_tcs_aoclc_simulator** to support the AOS closed loop simulation.
* Put the telescope related classes into the module of **telescope**.

.. _lsst.ts.phosim-1.1.2:

-------------
1.1.2
-------------

* Use the ``eups``, ``documenteer``, and **plantUML**.
* Use the **ts_wep** module.
* Use the scientific pipeline of ``sims_w_2019_18``.

.. _lsst.ts.phosim-1.1.1:

-------------
1.1.1
-------------

* Updated to use the scientific pipeline of ``sims_w_2019_02``.
* Reuse the **FilterType** Enum from **ts_tcs_wep**.

.. _lsst.ts.phosim-1.1.0:

-------------
1.1.0
-------------

* Refactor the code to decrease the number of function inputs.

.. _lsst.ts.phosim-1.0.0:

-------------
1.0.0
-------------

* Update the information and add the example scripts.

.. _lsst.ts.phosim-0.1.0:

-------------
0.1.0
-------------

* Initially integrate WEP and PhoSim.
