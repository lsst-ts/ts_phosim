.. py:currentmodule:: lsst.ts.phosim

.. _lsst.ts.phosim-modules:

##########
Modules
##########

The classes and files for each module are listed below.

.. _lsst.ts.phosim-modules_phosim:

-------------
phosim
-------------

This module is a high-level module to use other modules.

.. uml:: uml/phosimClass.uml
    :caption: Class diagram of phosim

* **PhosimCmpt**: High-level class to use the module of ts_phosim.
* **OpdMetrology**: OPD related metrology.
* **SkySim**: Sky simulator to add the stars.
* **PlotUtil**: Plot utility functions.
* **Utility**: Enums and functions used in this module.
* **MetroTool**: Metrology related functions contain the atmosphere model.
* **CloseLoopTask**: Close loop task to run the simulation with PhoSim.

.. _lsst.ts.phosim-modules_phosim_telescope:

-------------------
phosim.telescope
-------------------

.. uml:: uml/telescopeClass.uml
    :caption: Class diagram of phosim.telescope

* **PhosimCommu**: Interface to PhoSim.
* **CamSim**: Camera distortion correction.
* **MirrorSim**: Parent class of M1M3Sim and M2Sim classes.
* **M1M3Sim**: M1M3 mirror distortion of gravity and temperature gradient.
* **M2Sim**: M2 mirror distortion of gravity and temperature gradient.
* **TeleFacade**: Telescope facade pattern that intergate the correction of camera and mirror distortion correction to PhoSim.
