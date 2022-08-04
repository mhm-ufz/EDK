# -*- coding: utf-8 -*-
"""
Clean up the edk output in the check case folders.

Author
------
    Sebastian Mueller

Created
-------
    Jun 2022
"""
import os
import shutil
import glob

# output folder and reference folder
OUT = "output"


if __name__ == "__main__":
    # checking path
    cases_path = os.path.dirname(os.path.realpath(__file__))
    # get all cases folders (in the cases_path)
    cases = glob.glob(os.path.join(cases_path, "case*"))
    # sort the cases by name
    cases.sort()
    for case in cases:
        base = os.path.basename(case)
        print(base)
        # get the output directory
        out_dir = os.path.join(case, OUT)
        # remove it
        shutil.rmtree(out_dir, ignore_errors=True)
        logs = glob.glob(os.path.join(case, "*_" + base + "_log.txt"))
        for log in logs:
            os.remove(log)
    # remove common log dir
    shutil.rmtree("logs", ignore_errors=True)
    # final result
    print(" ..cleaned.")
