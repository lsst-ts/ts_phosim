# WEP Integrator with PhoSim

*This module is to integrate the WEP, PhoSim, and DAQ test stand. The original work was done by Bo Xian and Chuck Claver. The source code can be found in:*
[IM link](https://github.com/bxin/IM).

## 1. Version History

*Version 0.1*
<br/>
*Initially integrate WEP and PhoSim.*
<br/>
<br/>
*Version 1.0*
<br/>
*Update the information and add the example scripts.*
<br/>
<br/>
*Version 1.1.0*
<br/>
*Refactor the code to decrease the number of function inputs.*
<br/>
<br/>
*Version 1.1.1*
<br/>
*Updated to use the scientific pipeline of sims_w_2019_02. Reuse the FilterType Enum from ts_tcs_wep.*
<br/>

*Author: Te-Wei Tsai*
<br/>
*Date: 01-24-2019*

## 2. Platform

- *CentOS 7*
- *python: 3.6.6*
- *scientific pipeline (newinstall.sh from master branch)*
- *phosim_syseng4 (branch: aos, tag: firstdonuts)*

## 3. Needed Package

- *lsst_sims (-t sims_w_2019_02)*
- *ts_tcs_wep - develop branch (commit: b5dcb9a)*

## 4. Use of Module

*1. Setup the DM environment:*
<br/>
`source $path_of_lsst_scientific_pipeline/loadLSST.bash`
<br/>
`setup sims_catUtils -t $user_defined_tag -t sims_w_2019_02`
(e.g. `setup sims_catUtils -t ttsai -t sims_w_2019_02`)

*2. Setup the WEP environment:*
<br/>
`export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_wep/python`
<br/>
(e.g. `export PYTHONPATH=$PYTHONPATH:/home/ttsai/Document/stash/ts_tcs_wep/python`)

*3. Setup the wepPhoSim environment:*
<br/>
`export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_wep_phosim/python`
<br/>
(e.g. `export PYTHONPATH=$PYTHONPATH:/home/ttsai/Document/stash/ts_tcs_wep_phosim/python`)

## 5. Content
*This module contains the following classes and functions ([class diagram](./doc/wepPhosimClassDiag.png)):*

- **PhosimCommu**: Interface to PhoSim.
- **MetroTool**: Metrology related functions contain the atmosphere model.
- **OpdMetrology**: OPD related metrology.
- **CamSim**: Camera distortion correction.
- **MirrorSim**: Parent class of M1M3Sim and M2Sim classes.
- **M1M3Sim**: M1M3 mirror distortion of gravity and temperature gradient.
- **M2Sim**: M2 mirror distortion of gravity and temperature gradient.
- **TeleFacade**: Telescope facade pattern that intergate the correction of camera and mirror distortion correction to PhoSim.
- **SkySim**: Sky simulator to add the stars.
- **PlotUtil**: Plot utility functions.
- **Utility**: Enums and functions used in this module.

## 6. Example Script

- **calcOpd.py**: Test the OPD without the subsystem perturbation.
- **calcOpdAndSubSys.py**: Test the OPD with the subsystem perturbation.
- **checkStarAndSubSys.py**: Test the star donut in LSST camera with the subsystem perturbation.
- **checkWfsStarCoor.py**: Test to add the stars on WFS and get the images.
- **checkStarCoor.py**: Test to add the star by pixel position and get the image.
- **checkStarCoorWiLsstFAM.py**: Test to add the star by pixel position in LSST FAM condition and get the images.
