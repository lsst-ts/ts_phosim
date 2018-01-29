# WEP Integrator with PhoSim

*This module is to integrate the WEP, PhoSim, and DAQ test stand. The original work was done by Bo Xian and Chuck Claver. The source code can be found in:*

[IM link](https://github.com/bxin/IM)

## 1. Version History

*Version 0.1*
<br/>
*Initially integrate WEP and PhoSim.*

*Author: Te-Wei Tsai*
<br/>
*Date: 12-29-2017*

## 2. Platform

- *python: 3.6.2*
- *scientific pipeline: v14*
- *phosim_syseng2*

## 3. Needed Package

*N/A*

## 4. Use of Module

*N/A*

## 5. Content

*PhosimCommu.py: Interface to PhoSim.*
<br/>
*MetroTool.py: Metrology related functions contain the atmosphere model.*
<br/>
*OpdMetrology.py: OPD related metrology.*
<br/>
*CamSim.py: Camera distortion correction.*
<br/>
*MirrorSim.py: Parent class of M1M3Sim and M2Sim classes.*
<br/>
*M1M3Sim.py: M1M3 mirror distortion of gravity and temperature gradient.*
<br/>
*M2Sim.py: M2 mirror distortion of gravity and temperature gradient.*
<br/>
*CoTransform.py: Coordination transformation functions.*
<br/>
*TeleFacade.py: Telescope facade pattern that intergate the correction of camera and mirror distortion correction to PhoSim.*
<br/>
*SkySim.py: Sky simulator to add the stars.*
<br/>
*WEPController.py: Wavefront estimation controller class. This is a high level class to use the wavefront estimation pipeline.*

## 6. Example Script

*testOpd.py: Test the OPD without the subsystem perturbation.*
<br/>
*testOpdAndSubSys.py: Test the OPD with the subsystem perturbation.*
<br/>
*testQueryDbCoor.py: Test to add the stars by querying the UW BSC and get the images.*
<br/>
*testStarAndSubSys.py: Test the star donut in LSST camera with the subsystem perturbation.*
<br/>
*testStarAndSubSysWiComCam.py: Test the star donut in ComCam with the subsystem perturbation.*
<br/>
*testStarCoor.py: Test to add the star by pixel position and get the image.*
<br/>
*testWfsStarCoor.py: Test to add the stars on WFS and get the images.*