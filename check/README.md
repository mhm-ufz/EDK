# edk check cases

With the provided python script `edk_checks_run.py` you can run all checks in the `case_**` folders with one or more given edk executable(s).

It will be checked, if edk finishes normally and if the output matches the reference values.

A summary for all cases is given at the end.

## Usage

    python edk_checks_run.py [-h] [-e EXE [EXE ...]] [-v] [-l LOG_PATH] [-t OPENMP_THREADS] [-s [SKIP [SKIP ...]]] [-o [ONLY [ONLY ...]]]

Run the edk check cases with a given edk executable.

## Optional Arguments

    -h, --help            show this help message and exit
    -e EXE [EXE ...], --exe EXE [EXE ...]
                          Paths to edk exe[s]. (default: ['../edk'])
    -v, --verbose         Show the edk output. (default: no)
    -l LOG_PATH, --log_path LOG_PATH
                          Directory for edk-logs. (default: the resp. case dir)
    -t OPENMP_THREADS, --threads OPENMP_THREADS
                          Number of threads for openMP. No mpi allowed! (default: 0)
    -s [SKIP [SKIP ...]], --skip [SKIP [SKIP ...]]
                          skip cases (case_01 case_03 ..) (default: [])
    -o [ONLY [ONLY ...]], --only [ONLY [ONLY ...]]
                          only run cases (case_01 case_03 ..) (default: all)

## Examples

Run edk from parent directory in verbosity mode with mpi on 4 processes:

        python edk_checks_run.py -e ../edk -v -m 4

Silently run edk (given with an absolute path) with openmp on 4 threads:

        python edk_checks_run.py -e /abspath/edk_openmp -t 4

Silently run edk from parent directory:

        python edk_checks_run.py

Run with multiple edk exes:

        python edk_checks_run.py -e ../edk1 ../edk2

## Cleanup

To remove the created output of edk run:

        python edk_checks_clean.py

## Author

    Sebastian Mueller

## Contributors

    Stephan Thober, Robert Schweppe

Written Jan. 2022.
