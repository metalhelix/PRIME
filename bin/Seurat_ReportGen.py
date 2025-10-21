#!/usr/bin/env python

import os
import subprocess
import sys
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--rmd')
parser.add_argument('--cellbender_folder')
parser.add_argument('--mt_file')
args=parser.parse_args()

second_lib_cellbender_folder = args.cellbender_folder
rmd_file = args.rmd
mt_list_file = args.mt_file

shutil.copy2(rmd_file,second_lib_cellbender_folder)

os.chdir(second_lib_cellbender_folder)

result = subprocess.run(
    [f"R -e \"rmarkdown::render('Seurat_Report_gen.rmd', params=list(Reference_mito_list='{mt_list_file}'))\""],
    shell=True,
    capture_output=True,
    text=True,
)

if result.returncode == 0:
    print("Seurat report generated successfully!")
    subprocess.run(["mkdir", "PRIME_SC_out"])
    subprocess.run(
        [
            "mv",
            "Seurat_Report_gen.html",
            f"PRIME_SC_out/SeuratReport.html",
        ]
    )
    subprocess.run(
        [
            "mv",
            "SeuratObj.rds",
            f"PRIME_SC_out/SeuratObj.rds",
        ]
    )
    subprocess.run(
        [
            "mv",
            "Seurat_Report_gen_plots",
            f"PRIME_SC_out/Seurat_Report_gen_plots",
        ]
    )

    subprocess.run(["mv","-f" ,"PRIME_SC_out", ".."])
else:
    # If there was an error, print the error message
    sys.exit(f"Error: {result.stderr} ")
