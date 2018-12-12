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
*Version 1.1*
<br/>
*Refactor the code to decrease the number of function inputs.*

*Author: Te-Wei Tsai*
<br/>
*Date: 12-12-2018*

## 2. Platform

- *python: 3.6.6*
- *scientific pipeline: v16*
- *phosim_syseng2*

## 3. Needed Package

- *lsst_sims (-t sims_w_2018_47)*
- *ts_tcs_wep*

## 4. Use of Module

*1. Setup the DM environment:*
<br/>
`source $path_of_lsst_scientific_pipeline/loadLSST.bash`
<br/>
`setup sims_catUtils -t $user_defined_tag -t sims_w_2018_47`
(e.g. `setup sims_catUtils -t ttsai -t sims_w_2018_47`)

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

- **PhosimCommu**: Interface to PhoSim.
- **MetroTool**: Metrology related functions contain the atmosphere model.
- **OpdMetrology**: OPD related metrology.
- **CamSim**: Camera distortion correction.
- **MirrorSim**: Parent class of M1M3Sim and M2Sim classes.
- **M1M3Sim**: M1M3 mirror distortion of gravity and temperature gradient.
- **M2Sim**: M2 mirror distortion of gravity and temperature gradient.
- **CoTransform**: Coordination transformation functions.
- **TeleFacade**: Telescope facade pattern that intergate the correction of camera and mirror distortion correction to PhoSim.
- **SkySim**: Sky simulator to add the stars.
- **PlotUtil**: Plot utility functions.
- **Utility**: Enums and functions used in this module.

## 6. Example Script

- **testOpd.py**: Test the OPD without the subsystem perturbation.
- **testOpdAndSubSys.py**: Test the OPD with the subsystem perturbation.
- **testStarAndSubSysWiComCam.py**: Test the star donut in ComCam with the subsystem perturbation.
- **testStarAndSubSys.py**: Test the star donut in LSST camera with the subsystem perturbation.
- **testQueryDbCoor.py**: Test to add the stars by querying the UW BSC and get the images.
- **testWfsStarCoor.py**: Test to add the stars on WFS and get the images.
- **testStarCoor.py**: Test to add the star by pixel position and get the image.
- **testStarCoorWiComCam.py**: Test to add the star by pixel position in ComCam condition and get the images.
- **testStarCoorWiLsstFAM.py**: Test to add the star by pixel position in LSST FAM condition and get the images.
- **testWfsStarCoorAll.py**: Test to add the stars by pixel position for all corner WFS and get the images.