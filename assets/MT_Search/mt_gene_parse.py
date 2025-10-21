#!/usr/bin/env python 

import argparse
import pandas as pd
import os
import sys 

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_gtf')
parser.add_argument('-c', '--chrID')
parser.add_argument('-o', '--outfolder', default="./")
args=parser.parse_args()

# in_file = "/n/analysis/genomes/Artibeus_jamaicensis/GCF_021234435.1/annotation/NCBI_2023/gtfs/GCF_021234435.1.NCBI_2023.gtf"
# chrID = "NC_002009.1"

in_file = args.in_gtf
chrID = args.chrID

df_gtf = pd.read_csv(in_file, sep="\t", comment='#', header=None)
df_ = df_gtf[df_gtf[0] == chrID]

gene_name_feature_sec = "gene_id"
gene_name_feature_prime = "gene_name"

def gene_feature_fetch(input_value):
    feature_list = [v for v in input_value.split(";") if gene_name_feature_prime in v]
    
    if len(feature_list) == 0:
        feature_list = [v for v in input_value.split(";") if gene_name_feature_sec in v]
        feature_name = feature_list[0].replace(gene_name_feature_sec,"").strip().replace('"','')
    else:
        feature_name = feature_list[0].replace(gene_name_feature_prime,"").strip().replace('"','')

    print(feature_list)
    
    return feature_name

df_["featureName"] = df_[8].apply(gene_feature_fetch)

df_mt = pd.DataFrame()
df_mt[0] = df_["featureName"].unique()
df_mt.to_csv(f"{args.outfolder}/{os.path.basename(in_file).replace('.gtf','')}.mtGeneList.txt", index=False, header=False)