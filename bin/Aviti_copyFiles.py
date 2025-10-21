#!/usr/bin/env python

import glob
import argparse
import shutil
import os 

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--fastq_dir')
parser.add_argument('-l', '--libID')
parser.add_argument('-n', '--nanalysis_dir', help='the /n/analysis dir that the fastqs will be copied to')
args=parser.parse_args()


newR1_fastq_file = f"{args.fastq_dir}/{args.libID}_S*_L00*_R1_*.fastq.gz"
newR2_fastq_file = f"{args.fastq_dir}/{args.libID}_S*_L00*_R2_*.fastq.gz"


#Copy renamed fastqs to /n/analysis folders
if os.path.exists(args.nanalysis_dir) == False:
    os.makedirs(args.nanalysis_dir)

if os.path.isfile(os.path.join(args.nanalysis_dir, os.path.basename(newR1_fastq_file))) == False:
    shutil.copy2(newR1_fastq_file, args.nanalysis_dir)

if os.path.isfile(os.path.join(args.nanalysis_dir, os.path.basename(newR2_fastq_file))) == False:
    shutil.copy2(newR2_fastq_file, args.nanalysis_dir)
