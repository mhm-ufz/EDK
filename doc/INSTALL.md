# Compilation

The section 'Dependencies' lists the general requirements
for the compilation. The section 'System-dependent dependency installation'
gives some instructions on how to install these dependencies on Windows,
some Linux-distributions and MacOS.
Conda based dependency installation is described in the section 'Conda dependent installation'.


## Dependencies

To clone and compile EDK you need at least the following:

* Fortran compiler: We support [gfortran](https://gcc.gnu.org/fortran/), [nagfor](https://www.nag.com/content/nag-fortran-compiler) and [ifort](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
* Build system: We support [make](https://www.gnu.org/software/make/) and [ninja](https://ninja-build.org/)
* [cmake](https://cmake.org/): Software for build automation
* [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran): NetCDF I/O for Fortran
* [Lapack](https://www.netlib.org/lapack/): the linear algebra package
* [git](https://git-scm.com/): version control system
* (optional) [fypp](https://github.com/aradi/fypp): Fortran pre-processor written in Python


## System-dependent dependency installation

After you installed all dependencies on your system you can proceed with cloning and compiling.


### Unix (Linux / MacOS)

1. MacOS with [homebrew](https://brew.sh) available:
    ```bash
    brew install git gcc netcdf lapack cmake
    ```

2. Ubuntu, Mint and other apt-get based systems with matching repositories:
    ```bash
    sudo apt-get install git gfortran netcdf-bin libnetcdf-dev libnetcdff-dev liblapack-dev cmake
    ```

3. Archlinux:
    ```bash
    sudo pacman -S git gcc-libs netcdf-fortran lapack cmake
    ```

4. yum based systems (CentOS, OpenSuse):
    ```bash
    sudo yum -y install git gcc-gfortran netcdf-fortran lapack cmake
    ```


### Windows

On Windows 10 and later we recommend to use [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (WSL) to be able to use Linux by e.g. [installing Ubuntu](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10) there.

Easiest way to do so is:

1. install the [Windows Terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701)
2. open the Windows Terminal and type:
    ```bash
    wsl --install -d ubuntu
    ```

3. Open Ubuntu from the new entry in the start menu

Then you can follow the install instructions for Ubuntu from above.


### Module systems

If you are on a module system, load the modules gcc or intel depending on your
favorite compiler. Then, load the modules netcdf-fortran and cmake.

These modules will have system specific names, environments, etc.
You may use `module spider` to find the right packages and the
right dependencies, potentially use corresponding wiki pages.


#### On EVE (the cluster at the UFZ)

A set of load-scripts is provided in `hpc-module-loads` (see the [repository](https://git.ufz.de/chs/HPC-Fortran-module-loads) for more details), to load all needed modules for specific compilers:

- Example: GNU 7.3 compiler (`foss/2018b` Toolchain):
  ```bash
  source hpc-module-loads/eve.gcc73
  ```
  or (MPI support)
  ```bash
  source hpc-module-loads/eve.gcc73MPI
  ```


## Conda dependent installation

The simplest way to compile this project on your local computer is to use a [conda](https://docs.conda.io/en/latest/) environment (on Linux (including WSL) or MacOS) provided by [Miniforge](https://github.com/conda-forge/miniforge) to install NetCDF, a Fortran compiler, make, cmake.

On Windows we recommend to use [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (WSL) to be able to use Linux (see above) and set up conda there.

You can get the latest version of **Miniforge** with (for Linux/MacOS/WSL):
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
bash Miniforge3-$(uname)-$(uname -m).sh
```

To create a (local) conda environment with all dependencies type the following:
```bash
conda create -y --prefix ./fortran_env
conda activate ./fortran_env
conda install -y git cmake make fortran-compiler netcdf-fortran
```

Then you can proceed with cloning and compiling.


## Cloning the repository

First you need to clone the repository (if you already have `git`, otherwise see below):
```bash
git clone https://git.ufz.de/chs/progs/edk_nc.git
```

This will give you a new folder `edk_nc/` containing the whole repository. You can go into it by:
```bash
cd edk_nc/
```

If you then want to compile a specific version (different from the latest development version), you can check that out with e.g.:
```bash
git checkout v3.0.0
```

Afterwards you can continue with the compilation.


## Compilation commands

It could be necessary to set your desired fortran compiler with an environment variable, e.g.:
```bash
export FC=gfortran
```

We prepared a set of scripts to automatize the build and compilation process to generate an executable in the root directory with the following naming scheme:

- Release version `edk`:
  ```bash
  source scripts/compile
  ```
- Debug version `edk_debug`:
  ```bash
  source scripts/compile_debug
  ```
- Release version with OpenMP support `edk_openmp`:
  ```bash
  source scripts/compile_OpenMP
  ```
- Debug version with OpenMP support `edk_openmp_debug`:
  ```bash
  source scripts/compile_OpenMP_debug
  ```

Then you can find an executable `edk` (or `edk[_openmp][_debug]`) in the current folder.
You can execute it with:
```bash
./edk
```


## Installation

To install edk after compilation, i.e. make it available as a command `edk`, you can do the following (assuming you used the release compile script, otherwise replace `release` with the respective build folder):
```bash
cmake --install release
```

If you need to provide a prefix, where to install it, you can just pass one. For example, if you used a conda environment for compilation, you can also install edk there with:
```bash
cmake --install release --prefix $CONDA_PREFIX
```


## Compilation without Internet

Starting with version 3.0, EDK is depending on [FORCES](https://git.ufz.de/chs/forces/), our Fortran library for Computational Environmental Systems.
This library is downloaded on the fly by [CPM](https://github.com/cpm-cmake/CPM.cmake), the cmake package manager.

If you don't want to download it indirectly, know you wont have internet during your development or you want to work on routines provided by FORCES, you can place a copy of the FORCES repository in the root of your cloned EDK repository by e.g.:
```bash
git clone https://git.ufz.de/chs/forces.git
```
The new folder `forces/` will be automatically recognized during compilation as described above and nothing will be downloaded.

If you just want a specific version (see `src/CMakeLists.txt` for the currently used one), do this:
```bash
git clone --branch v0.3.2 --depth 1 https://git.ufz.de/chs/forces.git
```

If you have already cloned FORCES somewhere else, you can also provide a path to this repository. You can do this with all mentioned compile scripts, e.g.:
```bash
source scripts/compile -DCPM_forces_SOURCE=<path/to/your/forces/repo>
```

For example, if you have cloned FORCES next to `edk_nc`, this could look like this:
```bash
source scripts/compile -DCPM_forces_SOURCE=../forces
```


## Additional CMake infos

The presented compile scripts all just execute two cmake commands with a specific set of configuration flags.
The basic cmake workflow, to configure and compile in a `build/` folder, is:
```bash
cmake -B build
cmake --build build
```

You can control all `cmake` options by passing them as directives staring with `-D` to the cmake configuration.
For example for debug configuration, you can do the following:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
```

To configure the build interactively, you can also use [ccmake](https://cmake.org/cmake/help/latest/manual/ccmake.1.html) (command line tool) or the [CMake GUI](https://cmake.org/cmake/help/latest/manual/cmake-gui.1.html) (graphical user interface).
Check their respective documentation for further information.
To use `ccmake` you can do the following:
```bash
cmake -B build
ccmake build
```

Then set your desired options and re-configure your build (by pressing `c`).
Afterwards build you project as always by executing:
```bash
cmake --build build --parallel
```
