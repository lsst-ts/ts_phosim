.. py:currentmodule:: lsst.ts.phosim

.. _lsst.ts.phosim-version_history:

##################
Version History
##################

.. _lsst.ts.phosim-1.1.7:

-------------
1.1.7
-------------

Use sims_w_2019_31. Use the latest versions of ts_wep and ts_ofc. Remove the conda package installation in Jenkinsfile. Update the permission of workspace after the unit test.

.. _lsst.ts.phosim-1.1.6:

-------------
1.1.6
-------------

Use sims_w_2019_29. Supress the warning in unit tests. Fix the warning of nan in atmosphere structure function. Rotate the OPD and support the sky file, minimum DOF, and M1M3 force error ratio in command line tasks.

.. _lsst.ts.phosim-1.1.5:

-------------
1.1.5
-------------

Use sims_w_2019_24. Support the eimage in comcamCloseLoop.py. Depend on the SensorWavefrontError in ts_wep. Update the table file.

.. _lsst.ts.phosim-1.1.4:

-------------
1.1.4
-------------

Minor bugs fixed. Add the get methods for SkySim and OpdMetrology classes. Use the CamType of ts_wep module in TeleFacade class. Update PhosimCmpt class to use the interface classes of ts_wep and ts_ofc. Use the scientific pipeline of sims_w_2019_20. Add the command line tasks of close-loop simulation.

.. _lsst.ts.phosim-1.1.3:

-------------
1.1.3
-------------

Combine with ts_tcs_aoclc_simulator to support the AOS closed loop simulation. Put the telescope related classes into the module of telescope.

.. _lsst.ts.phosim-1.1.2:

-------------
1.1.2
-------------

Use the eups, documenteer, and plantUML. Use the ts_wep module. Use the scientific pipeline of sims_w_2019_18.

.. _lsst.ts.phosim-1.1.1:

-------------
1.1.1
-------------

Updated to use the scientific pipeline of sims_w_2019_02. Reuse the FilterType Enum from ts_tcs_wep.

.. _lsst.ts.phosim-1.1.0:

-------------
1.1.0
-------------

Refactor the code to decrease the number of function inputs.

.. _lsst.ts.phosim-1.0.0:

-------------
1.0.0
-------------

Update the information and add the example scripts.

.. _lsst.ts.phosim-0.1.0:

-------------
0.1.0
-------------

Initially integrate WEP and PhoSim.
