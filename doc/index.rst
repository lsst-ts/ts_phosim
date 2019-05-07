.. py:currentmodule:: lsst.ts.phosim

.. _lsst.ts.phosim:

##############
lsst.ts.phosim
##############

This module is a high-level module to use PhoSim. This can put the perturbations of camera temperature, mirror gravity, and mirror temperature corrections into PhoSim. This can also put the perturbation of telescope degree of freedom (DOF) into PhoSim. This can combine with WEP and OFC module to do the closed-loop simulation on PhoSim as a test framework.

.. _lsst.ts.phosim-using:

Using lsst.ts.phosim
====================

.. toctree::
   :maxdepth: 1

Important classes:

* `TeleFacade` intergates the correction of camera and mirror distortion correction to PhoSim.
* `OpdMetrology` does the optical path difference (OPD) related metrology.
* `SkySim` adds the stars.

.. _lsst.ts.phosim-pyapi:

Python API reference
====================

.. automodapi:: lsst.ts.phosim
    :no-inheritance-diagram:

.. _lsst.ts.phosim-content:

Content
====================

.. toctree::

   content

.. _lsst.ts.phosim-contributing:

Contributing
============

``lsst.ts.phosim`` is developed at https://github.com/lsst-ts/ts_phosim.

.. _lsst.ts.phosim-version:

Version
====================

.. toctree::

   versionHistory
