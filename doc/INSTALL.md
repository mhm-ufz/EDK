Dependencies:
=============
For Windows, some Linux distributions and soon also MacOS
specific installation instructions for the following list
can be found below.

- a fortran compiler
- make (a tool to compile a program)
- cmake (version >= 3.5) (a tool to create a system dependent makefile)
- fitting netcdf-fortran libraries (libraries for the usage of the data format netcdf on which edk depends)
- (optional, but makes things much easier) git

Git is a version-control system. If you want to contribute to a project, it is highly recommended to
use Git. You can use Git to download (clone) the project to your local pc and have a look at the history or
synchronize it without copying the whole repository again. You can also download the project folder without
Git, but this would not allow you to pull updates from and push changes to our repository.

System dependend installation instructions:
===========================================
### Windows:
[Cygwin](https://cygwin.com/) is an environment with a terminal that allows to compile and
run programs of Unix-like systems. You can find further instructions to install cygwin on the webpage, as well as
instructions on how to install further dependencies after the installation.

After the installation of cygwin and its dependencies edk will be installed
using cygwin. All commands and the execution of edk only run in that environment.

Install cygwin by executing the cygwin setup and choose the following dependencies:

- [ ] gcc-fortran (the fortran compiler)
- [ ] make
- [ ] cmake (version >= 3.5)
- [ ] libnetcdf-fortran-devel
- [ ] Git *(optional, Git is also available outside of cygwin, [see the Git website](https://git-scm.com/downloads))*

While installing cygwin you will have to choose a mirror. A mirror is a server
on the internet where the files for the installation come from. Choose any server
located near your city and when in doubt, choose the first one in the list.
In the next step you can find all available packages provided by cygwin, set
the view to "full". In the search panel you can filter the packages 
by the dependencies listed above (e.g. make). When you choose a
version, the newest one is usually a good choice if not marked as experimental.

*Note for UFZ members:* Install cygwin locally, do not choose a location on the
network for the installation.

Some cygwin versions create a new home directory for you. You may check e.g. here:

    C:\cygwin64\home\$username


### Ubuntu, Mint and other apt-get based systems with matching repositories:

    sudo apt-get install git # (optional)
    sudo apt-get install gfortran netcdf-bin libnetcdf-dev libnetcdff-dev cmake libblas-dev liblapack-dev

### Archlinux:

    sudo pacman -S git # (optional)
    sudo pacman -S gcc-libs netcdf-fortran cmake

### Module systems:

If you are on a module system, load the modules gcc or intel depending on your
favorite compiler. Then, load the modules netcdf-fortran and cmake. 

These modules will have system specific names, environments, etc.
You may use `module spider` to find the right packages and the
right dependencies, potentially use corresponding wiki pages.

#### On eve (the cluster at the UFZ):

From 2020 onward the software is installed in toolchains, where all included software specific compilers
From the source directory use a script provided in `moduleLoadScripts`,
for example for the GNU 7.3 compiler:

    source moduleLoadScripts/eve.gfortran73
    
for intel use
    source moduleLoadScripts/eve.intel[version]
    
    ATTENTION!: in intel cmake cannot find LAPACK and BLAS Libraries, the paths need to be set manually by
    cmake -DLAPACK_LIBRARIES=$MKLROOT/lib/intel64/libmkl.so -DBLAS_LIBRARIES=$MKLROOT/lib/intel64/libmkl.so ..

### MacOS:

*(to be added)*

Specific setups:
================

The following hints can replace the step `cmake ..` in the installation instruction.

You can skip this part and continue with "Installation", if you do not have a module system
setup (like on clusters) or if you have not installed all packages with a package manager,
such as cygwin or apt-get.

### Module systems:

The executable can be build in a way that it runs independend of loaded modules in the end. The
module system, though, adds system paths in the backround the user should not care about too much, so
the setup is a workaround. (This would be the case with any other building tool aswell.)
It should be stable, anyway.

In case you want to have a module-independend build, instead of just executing `cmake ..`, either run

    cmake -DCMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND:STRING=ON ..

or

    cmake -C ../CMakeCacheFiles/eve ..

or change the variable `CMAKE_BUILD_MODULE_SYSTEM_INDEPENDEND` with `ccmake` to `ON` after running `cmake ..`.

### None standard locations for the netcdf-library (e.g. standard setup Macs in CHS):

Find the location of the `nf-config` file, for example by using:

    find / -iname "*nf-config*" 2>/dev/null

This searches the root directory `/` for a file with a name containing the string "nf-config", not
taking into account upper and lower case. It writes error messages like "permission denied" into
the void.

Then, instead of running `cmake ..` if not using the standard compiler,
set the fortran compiler variable to the wished compiler, e.g.

    export FC=gfortran

then either run

    cmake -DCMAKE_NETCDF_DIR:STRING=/path/to/nf-config/of/used/compiler

or copy the file `specificSetup` to some other place:

    cp ../CMakeCacheFiles/specificSetup .

However, in case you want to keep it, you should choose a place
outside the build repository. Edit the file as follows:
add the path to your `nf-config` file, and after editing, run:

    cmake -C specificSetup ..

or change the variable `CMAKE_NETCDF_DIR` to the path to the `nf-config` file with `ccmake` after running `cmake ..`.

Installation
============

1. Change to a directory where you want to store the source code.
2. Clone the corresponding edk repository into a folder, either using Git (if installed):

        git clone -b cmake git@git.ufz.de:chs/progs/EDK.git edk/

3. Create and change to a build directory where you want to store the build, e.g. inside the Git source directory

        cd edk
        mkdir build

    Change into the build directory:

        cd build

4. Generate a system dependent makefile

    Execute `cmake` with the path to the Git source directory as parameter.

       cmake ..

    If everything worked well a Makefile was created with the corresponding paths.

    *Note: have a look at "Specific setups" above in case you are using module systems,
    or when the netcdf libraries are not located where the package manager usually installs libraries, 
    or when they are not saved in environment variables (i.e., classical MacOS setups at CHS).*

5. Make the build:

   Execute make:

        make

    If this also worked fine, an executable was created, which has to be moved or copied to the Git source directory.

6. Execute the file:

        cd ..
        cp build/edk .

    On Windows the executable is called `edk.exe` instead of `edk`. In that case
    instead of `cp build/edk .` execute 

        cp build/edk.exe .

    Now you might execute edk:

        ./edk

*Note concerning the development of the cmake setup: one could automatically
 link the executable with the `cmake` code inside the Git source directory
  which is not done for two reasons:*

- *The executable depends on your local system, so it should never be commited and pushed to other users.
    Nothing should be build inside the source directory which we did not do by hand.*
- *The directory where edk is executed usually is not the source directory but the directory where you want to run
   your tests. In case of the test setup it is the same, usually it is not.*

Building Realease or Debug versions:
====================================
If you want to set up specific versions of the build, you can
create different folders for that. Assume a release and a debug
version. Then a good idea would be to create one folder named
 `debug` and one folder named `release`

    mkdir release

    mkdir debug

inside the `release` folder one would execute

    cmake -DCMAKE_BUILD_TYPE=Release ..

and inside the `debug` folder

    cmake -DCMAKE_BUILD_TYPE=Debug ..

Executing

    make

in the corresponding folder would then always result in a release build or respectively in a debug build.

Trouble shooting:
=================

On brew/homebrew setup MacOS systems there is no working `nf-config` by now. Execute:

    nf-config --all

and if it says something like "is not implemented yet" the issue is not solved yet. But it is on my tracklist.

In any other case feel free to write an email to <mailto:maren.kaluza@ufz.de>.

**cmake** is far from being my main task, so it will probably take a while until I can track a problem.
I would be happy having bug reports, anyhow.
