.. py:currentmodule:: lsst.ts.phosim

.. _lsst.ts.phosim-phosimDonutTemplates:

########################################
Creating Phosim Donut Templates
########################################

This document describes the code to generate the phosim donut templates
that are used by `DonutTemplatePhosim` in `ts_wep`.

Running the code
================

1) Make sure your environment is set up to use ts_phosim and that you have phosim
   installed with the environment variable $PHOSIMPATH set to the phosim directory.
2) Use python to run `runCreatePhosimDonutTemplates.py` in the bin.src directory.
   Use `--help` to see all options.
3) Templates will appear in
   `ts_phosim/policy/donutTemplateData/phosimTemplates/`.

Input files
===========

The following input files found in `policy/donutTemplateData` are used
to create the donut templates.

* **starExtra.inst**: Phosim instance catalog to create one extra-focal donut per CCD
                      of LSST camera using PhosimMapper.
* **starIntra.inst**: Phosim instance catalog to create one intra-focal donut per CCD
                      of LSST camera using PhosimMapper.
* **star.cmd**: Phosim command file.
* **createPhosimDonutTemplateConfig.yaml**: Gen 3 Butler config file to run the ISR.
