# README for Daniele's subproject
## snakemake pipeline for limix GWAS on atmosphere

Using a medium2 Ubuntu 18.04 NoDesktop Base v4.0 atmosphere image with snakemake and conda installed
snakemake --use-conda --cores 1

## Step 1:
sudo -i or chown
cd /scratch


## Step 2: clone repository from github
git clone https://github.com/redtreevole/Foss2020Team6project

## Step 3: stage data in from Data Store at iplant/home/filiaultd/GWAS_data to directory 001.data
mkdir 001.data
cd 001.data
irsync -r i:GWAS_data ./

## Step 3: run snakemake pipeline
cd ..
snakemake --use-conda --cores 1
### I can't make the conda environment!!!

## Step 4: stageout results to Data Store
irsync -r /scratch/Foss2020Team6project/limixGWAS/001.data i:iplant/home/filiaultd/GWAS_data/003.results

## Other general notes
use cyberduck to move data to/from local computer to Data Store
