# PRIME

General PRIME documentation for development. 

---

## Running the pipeline

Refer to this page for general guidance on how to run the pipeline: [How to run PRIME](https://stowersinstitute.sharepoint.com/sites/ComputationalBiologyPublic/SitePages/PRIME.aspx)

---

## Code Structure

The source code for PRIME is stored under: /n/ngs/tools/PRIME/RUN/

In this and most nextflow pipelines:
- main.nf defines the general logic
- workflows/ contain all workflows defined in main.nf
- modules/ contain all processes defined in main and workflows under workflows/
- bin/ contain all the python and r scripts used in processes under module/ and workflows/
- scripts/ contain python scripts for pipeline kickoff (CronJobs) 
- nextflow.config contains all process parameters for slurm resource allocations 

---

## General Logic 

The parameters and the general logic of the pipeline is defined in "main.nf". 
In this file, I specified in a series of if-else statements on how to handel different order types and how to correspond workflows to order types.

***Add additional if-else statements in main.nf to include new order types further down the line.***

---

## Workflows

Workflows of different order types are designed to be independent of each other.
I would encourage people to add new workflows when adding new order types. 

---

## Modules

Processes are grouped based on order types. 
Processes can be shared between different workflows. And they are in PRIME to reduce redundancy. 
When adding a new workflow or order type, I would still recommend writing new processes. 

---

## CronJobs

CronJobs are used to automatically kick off prime.
CronJobs are ran under the compbio_svc account. 
Command: 
49 * * * * ~/miniconda3/envs/R-SECUNDO3/bin/python /n/ngs/tools/PRIME/RUN/scripts/CronJobRunner_PRIME_v4.py --samplesheet /n/analysis/genomes/sampleSheet_ROBOINDEX_2023.csv >> /n/ngs/tools/PRIME/logs/CronJob.log 2>&1

--- 

## Run Environment
PRIME is running within a conda environment stored here: /home/compbio_svc/miniconda3/envs/PrimeG4 

---

## Log Files

Log files for PRIME are stored under /n/ngs/tools/PRIME/logs
- CronJob.log : Direct output of the CronJob script
- PRIME_Orders.log | PACBIO_Orders.log : Orders that were detected in previous runs with time stamps
- machine_run_logs/ : Outputs of the pipeline per flowcell

---


