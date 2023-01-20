[![docs](https://img.shields.io/badge/docs-ts--phosim.lsst.io-brightgreen)](https://ts-phosim.lsst.io/)

# PhoSim Use

This module is a high-level module to use PhoSim.

See the docs: <https://ts-phosim.lsst.io/>

## Platform

- CentOS 7
- python: 3.7.8
- [phosim_syseng4](https://github.com/lsst-ts/phosim_syseng4) (branch: aos, tag: firstdonuts)

## Needed Package

- [ts_wep](https://github.com/lsst-ts/ts_wep)
- [ts_ofc](https://github.com/lsst-ts/ts_ofc)
- [phosim_utils](https://github.com/lsst-dm/phosim_utils)
- [black](https://github.com/psf/black) (optional)
- [documenteer](https://github.com/lsst-sqre/documenteer) (optional)
- [plantuml](http://plantuml.com) (optional)
- [sphinxcontrib-plantuml](https://pypi.org/project/sphinxcontrib-plantuml/) (optional)

## Use of Module

1. Setup the WEP and OFC environments first, and then setup and build phosim_utils:

```bash
cd $phosim_utils_directory
setup -k -r .
scons
```

2. Setup the PhoSim environment by `eups`:

```bash
cd $ts_phosim_directory
setup -k -r .
scons
```

3. Set the path variables:

```bash
export PHOSIMPATH=$path_to_phosim
export AOCLCOUTPUTPATH=$path_to_output
```

## Code Format

This code is automatically formatted by `black` using a git pre-commit hook.
To enable this:

1. Install the `black` Python package.
2. Run `git config core.hooksPath .githooks` once in this repository.

## Example Script

- **calcOpd.py**: Test the OPD without the subsystem perturbation.
- **calcOpdAndSubSys.py**: Test the OPD with the subsystem perturbation.
- **checkStarAndSubSys.py**: Test the star donut in LSST camera with the subsystem perturbation.
- **checkWfsStarCoor.py**: Test to add the stars on WFS and get the images.
- **checkStarCoor.py**: Test to add the star by pixel position and get the image.
- **checkStarCoorWiLsstFAM.py**: Test to add the star by pixel position in LSST FAM condition and get the images.

## Command Line Task

- **opdCloseLoop.py**: Close-loop simulation in the optical path difference (OPD) level, which means the wavefront estimation pipeline (WEP) is not considered.
- **imgCloseLoop.py**: Close-loop simulation of the images. This task supports the amplifier images and eimages of PhoSim.
- **createSkyFile.py**: Code to create input sky catalog files for the closed loop from reference catalogs available in a butler repository.

## Example Sky Files

There are two sky files in `tests/testData/sky/` directory that can be used in the test of command line task. One is for ComCam and one is for LSST full array mode (FAM): **skyComCam.txt** and **skyLsstFam.txt**.

## Build the Document

To build project documentation, run `package-docs build` to build the documentation.
The packages of **documenteer**, **plantuml**, and **sphinxcontrib-plantuml** are needed.
The path of `plantuml.jar` in `doc/conf.py` needs to be updated to the correct path.
To clean the built documents, use `package-docs clean`.
See [Building single-package documentation locally](https://developer.lsst.io/stack/building-single-package-docs.html) for further details.

## Reference of PhoSim with active optics (AOS)

- The original work was done by Bo Xin and Chuck Claver. The source code can be found in: [IM](https://github.com/bxin/IM).
