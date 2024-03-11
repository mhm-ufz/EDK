# The External Drift Kriging - EDK program

<div align="center">
<img src="https://git.ufz.de/chs/logos/-/raw/master/EDK.png" alt="EDK-LOGO" width="251px" style="width:251px;"/>
</div>

This repository contains the external drift kriging (EDK) Fortran program developed at the Dept. Computational Hydrosystems at the Helmholtz Centre for Environmental Research - UFZ.

The EDK program comes with a LICENSE agreement, this includes also the GNU Lesser General Public License.

**Please note**: The GitLab repository grants read access to the code.
If you like to contribute to the code, please contact stephan.thober@ufz.de or sebastian.mueller@ufz.de.

[TOC]

## Installation

Installation instructions can be found in [INSTALL](doc/INSTALL.md) for Windows, MacOS, and GNU/Linux distributions.

The simplest way to compile the EDK program is to use a [conda](https://docs.conda.io/en/latest/) environment (on Linux (including Windows/WSL) or MacOS)
provided by [Miniforge](https://github.com/conda-forge/miniforge):
```bash
conda create -y --prefix ./fortran_env
conda activate ./fortran_env
conda install -y git cmake make fortran-compiler netcdf-fortran liblapack
source scripts/compile
```
This will give an executable `edk`.

## Usage

To run the EDK program, you need a set of station data files, a look-up-table for these stations and a DEM file for the external drift.
All configuration is done with a namelist file, which is called `edk.nml` by default. See the example file for all input specifications.

Then you can just execute `./edk` next to this file.

You can also explicitly specify the namelist file and/or change the working directory by passing options the the `edk` command like:
```bash
./edk -n edk.nml check/case_01
```

To see the help text, execute:
```bash
./edk --help
```

## Cite as

Please refer to the EDK algorithm by citing Samaniego et al. (2011). EDK aplications in Samaniego et al. (2013) or Zink et al. (2017).
Rainfall network design and EDK cross-validation in Zacharias, S. et al. (2011).

- Samaniego, L., R. Kumar, and C. Jackisch (2011), "Predictions in a data-sparse region using a regionalized grid-based hydrologic model driven by remotely sensed data", Hydrology research, 42(5), 338–355, doi:10.2166/nh.2011.156.
- Samaniego, L. R. Kumar, M. Zink (2013), "Implications of Parameter Uncertainty on Soil Moisture Drought Analysis in Germany", J Hydrometeor, 2013 vol. 14 (1) pp. 47-68. http://journals.ametsoc.org/doi/abs/10.1175/JHM-D-12-075.1
- Zink, M., R. Kumar, M. Cuntz, and L. Samaniego (2017), "A high-resolution dataset of water fluxes and states for Germany accounting for parametric uncertainty", Hydrol. Earth Syst. Sci., 21(3), 1769–1790, doi:10.5194/hess-21-1769-2017.
- Zacharias, S., H. Bogena, L. Samaniego et al. (2011), "A Network of Terrestrial Environmental Observatories in Germany", Vadose Zone Journal, 10(3), 955, doi:10.2136/vzj2010.0139.

## License

LGPLv3 (c) 2005-2024 CHS-Developers
