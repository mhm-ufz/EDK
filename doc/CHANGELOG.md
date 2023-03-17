# Changelog

[TOC]

All notable changes to will be documented in this file.


## v3.0.0 - ?
See the git [diff](https://git.ufz.de/chs/progs/edk_nc/-/compare/v2.0.0...v3.0.0) for details.

### Enhancements
- documentation page
- command line interface:
  - set namelist file name and working directory: `edk -n edk.nml path/to/cwd`
  - get help message (`edk -h`) and version info (`edk -v`)
- added folder with compile scripts: `scripts/`
- new checking scripts for check cases (run in CI)
- added second check case with variogram fitting
- added version files `version.txt` and `version_date.txt`
- memory optimization ([5](https://git.ufz.de/chs/progs/edk_nc/-/merge_requests/5)):
  - the following global variables were replaced with a new class `edk_dist` that calculates these distances on demand to prevent memory issues:
    - `dCS` distance matrix between cells (C) and stations (S)
    - `dS` distance matrix between stations (S)
    - `dz2S` pairwise squared z-value differences betweens stations
- added option to invert y dimension (currently: invert_y = .True. in edk.nml for mHM compatible input)
- added feature to allow time buffering (setable in namelist)

### Changes
- use [FORCES](https://git.ufz.de/chs/forces/) as dependency
- move lonely routines to the modules where they are used
- renamed modules to unique naming scheme: `mo_edk_*`
- added variables `projection_name`, `variable_unit` and `variable_long_name` to netcdf output
- changed output variable to single precision

### Bugfixes
- function `tVar` now returns `real(dp)` instead of `real(8)`
- correction of wrong date allocation
