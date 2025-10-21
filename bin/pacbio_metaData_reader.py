#!/usr/bin/env python

import os
import sys
import argparse
import xml.etree.ElementTree as ET
import pandas as pd
import requests

requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS = 'ALL:@SECLEVEL=1'

def fetch_sample_table_from_LIMS(fc_id):
    """
    Get the run info from LIMS. Returns the /n/analysis/ folder where the primary analsysis pipeline copied all the data
    """
    # Get run info from lims
    NGS_LIMS = 'https://lims.stowers.org/zanmodules/molecular-biology/ngs'
    API_TOKEN = 'ca7952666a03dd4e59d0cd59e39fecc7' # Should get a new one for each pipeline, or even each user

    header = {'x-zan-apitoken': f'{API_TOKEN}', 
        'Accept': 'application/json'}
    run_info = requests.get(f'{NGS_LIMS}/flowcells/{fc_id}/samples', headers=header, verify=False)
    if not run_info.ok: 
        sys.exit('Malformed API request. Please double check your flowcell ID')
    lims_data = run_info.json() # Close request
    run_info.close() 
    df_run_info = pd.DataFrame.from_dict(lims_data)
    df_sample_info = pd.concat([df_run_info['readLength'] ,df_run_info['readType'], df_run_info['samples'].apply(pd.Series)], axis=1)
    return df_sample_info

def bamFinder(xml_path, barcode):
    fcpath = os.path.dirname(os.path.dirname(xml_path))
    for root, dirs, files in os.walk(os.path.join(fcpath, "hifi_reads")):
        for file in files:
            if barcode in file and file.endswith(".bam"):
                return os.path.join(root, file) if os.path.exists(os.path.join(root, file)) else None
    return None

def qcFinder(xml_path):
    fcpath = os.path.dirname(os.path.dirname(xml_path))
    for root, dirs, files in os.walk(os.path.join(fcpath, "statistics")):
        for file in files:
            if file.endswith(".report.pdf"):
                return os.path.join(root, file) if os.path.exists(os.path.join(root, file)) else sys.exit(f"Error: No qc report found in {os.path.join(fcpath, 'statistics')}")

parser = argparse.ArgumentParser(description="Extract metadata from PacBio Sequel II/IIe/Revio XML file and LIMS")
parser.add_argument("-x", "--xml", required=True, help="Path to the PacBio metadata XML file")
args = parser.parse_args()

xml_file = args.xml
qc_report = qcFinder(xml_file)
fc_path = os.path.dirname(os.path.dirname(xml_file))

xml_tree = ET.parse(xml_file)
root = xml_tree.getroot()
collections = root.findall(".//{*}Collections")[0] # There should be only one Collections tag in the xml file, i think.....

fcID = collections.find(".//{*}CellPac").attrib.get("LabelNumber")

"""Long Plex currently not supported in this pipeline"""
if "LongPlex" in collections.find(".//{*}WellSample").attrib.get("Description"):
    sys.exit("Error: LongPlex orders are not supported in this pipeline. Please contact the bioinformatics team for assistance.")
molng = collections.find(".//{*}WellSample").attrib.get("Description").split()[0]
""""""

libID_list = []
barcode_list = []

for item in collections.findall(".//{*}BioSample"): 
    libID = item.attrib.get("Name")
    barcodes = []
    for barcode in item.findall(".//{*}DNABarcode"): # There can only be 1 barcode per library
        barcode_name = barcode.attrib.get("Name")
        barcodes.append(barcode_name)
    if len(barcodes) > 1:
        sys.exit(f"Error: More than one barcode found for {libID}, {molng}")
    # print(f"Molng: {molng}, LibID: {libID}, Barcode: {barcodes[0]}")
    libID_list.append(libID)
    barcode_list.append(barcodes[0].split("--")[0])

df_ = pd.DataFrame({"Molng": [molng] * len(libID_list), "LibID": libID_list, "Barcode": barcode_list})
df_metadata = df_.drop_duplicates()
df_samplesheet = fetch_sample_table_from_LIMS(fcID)

df = pd.merge(df_metadata, 
                df_samplesheet[["libID", "resultPaths", "sampleName", "genomeVersion", "speciesName", "requestingDepartment", "orderType"]],
                left_on="LibID", right_on="libID", how="left").drop(columns=["libID"])

df["speciesName"] = df["speciesName"].apply(lambda s : s.split(",")[0].rstrip() if pd.notna(s) else s)
df["BamPath"] = df["Barcode"].apply(lambda x: bamFinder(xml_file, x) if x else sys.exit(f"Error: No bam file found for {x} in order: {molng}"))
df["QcReport"] = qc_report
df["FcPath"] = fc_path

df.to_csv("pacbio_metadata.csv", index=False)