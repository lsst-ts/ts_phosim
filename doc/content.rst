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

.. uml:: uml/phosimClass.uml
    :caption: Class diagram of phosim

* **PhosimCommu**: Interface to PhoSim.
* **MetroTool**: Metrology related functions contain the atmosphere model.
* **OpdMetrology**: OPD related metrology.
* **CamSim**: Camera distortion correction.
* **MirrorSim**: Parent class of M1M3Sim and M2Sim classes.
* **M1M3Sim**: M1M3 mirror distortion of gravity and temperature gradient.
* **M2Sim**: M2 mirror distortion of gravity and temperature gradient.
* **TeleFacade**: Telescope facade pattern that intergate the correction of camera and mirror distortion correction to PhoSim.
* **SkySim**: Sky simulator to add the stars.
* **PlotUtil**: Plot utility functions.
* **Utility**: Enums and functions used in this module.
