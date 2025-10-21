#!/usr/bin/env python

import argparse
import os
import pandas as pd
import json
import subprocess
import sys

parser = argparse.ArgumentParser()
# parser.add_argument("--molngID", "-m")
parser.add_argument("--secondaryPath", "-s")
args = parser.parse_args()

secondary_Path = args.secondaryPath
df_summary = pd.read_csv(os.path.join(secondary_Path, "scripts/summary.csv"))

df_summary["Seurat Obj"] = [
    f"{secondary_Path}/{v}/PRIME_SC_out/SeuratObj.rds"
    for v in df_summary["libID"].tolist()
]

lab = secondary_Path.split("secondary/")[-1].split("/MOLNG")[0].lstrip("/").split("/")[0]
requester = secondary_Path.split("secondary/")[-1].split("/MOLNG")[0].split("/")[1]
df_summary["lab"] = len(df_summary) * [lab]
df_summary["requester"] = len(df_summary) * [requester]

dic_summary = df_summary.to_dict(orient="records")

molngID_info = os.path.basename(os.path.dirname(f"{secondary_Path}/"))
dic_all = {}
dic_all.setdefault(molngID_info, dic_summary)
dic_json = {}
dic_json.setdefault(lab, dic_all)

json_file_name = f"{secondary_Path}/{molngID_info}.PRIMESC-DB.json"
with open(json_file_name, "w") as json_file:
    json.dump(dic_json, json_file, indent=4)

# Connect to BioTools

def run_cmd(command):
    try:
        # Run the subprocess and capture the output
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        # Print the output if successful
        print(f"Command {' '.join(command)} succeeded with output:\n{result.stdout}")
    except subprocess.CalledProcessError as e:
        # Handle errors (e.g., command failed)
        print(f"Command {' '.join(command)} failed with error:\n{e.stderr}")
    except Exception as e:
        # Handle any other unexpected exceptions
        print(f"An unexpected error occurred while running {' '.join(command)}:\n{e}")

run_cmd(["ssh", "compbiotools", "/opt/bioTools/db_script/venv/bin/python", "/opt/bioTools/db_script/insert_data.py", f"{os.path.join(os.getcwd(),json_file_name)}"])
