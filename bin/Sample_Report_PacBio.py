#!/usr/bin/env python

import pandas as pd
import argparse
import os
import sys

parser = argparse.ArgumentParser(description="Generate PacBio Sample_Report.csv")
parser.add_argument("-i", "--input", required=True, help="Input CSV file with PacBio metadata")
args = parser.parse_args()

df_samplesheet = pd.read_csv(args.input)

df = pd.DataFrame()

# e.g sample report csv:
# Output,Order,OrderType,SampleName,LibraryID,IndexSequence1,IndexSequence2,Read,Reference,Lab,TotalReads,AlignPercent,Type,ReadLength,Species,FastqPath
# m84222_241127_200255_s4.hifi_reads.bc2004.bam,MOLNG-4146,PacBio,S72385,L72385,,,,,,,,,,,/n/analysis/Gerton/md2370/MOLNG-4146-4160/EA156125/hifi_reads/m84222_241127_200255_s4.hifi_reads.bc2004.bam

df["Output"] = df_samplesheet["LibID"] + ".bam"
df["Order"] = df_samplesheet["Molng"]
df["OrderType"] = "PacBio"
df["SampleName"] = df_samplesheet["sampleName"]
df["LibraryID"] = df_samplesheet["LibID"]
df["IndexSequence1"] = df_samplesheet["Barcode"]
df["IndexSequence2"] = "NA"
df["Read"] = "NA"
df["Reference"] = df_samplesheet["genomeVersion"]
df["Lab"] = df_samplesheet["requestingDepartment"]
df["TotalReads"] = "NA"
df["AlignPercent"] = "NA"
df["Type"] = df_samplesheet["orderType"]
df["ReadLength"] = "NA"
df["Species"] = df_samplesheet["speciesName"].apply(lambda s : s.split("|")[0].rstrip() if pd.notna(s) else s)
df["FastqPath"] = df_samplesheet["resultPaths"] + df_samplesheet["LibID"] + ".bam"

df.to_csv("Sample_Report.csv", index=False)