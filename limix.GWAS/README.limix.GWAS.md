# README for Daniele's subproject
## snakemake pipeline for limix GWAS on atmosphere :snake:

Using a medium2 Ubuntu 18.04 NoDesktop Base v4.0 atmosphere image with snakemake and miniconda installed

## Step 1:
    sudo -i or chown
    cd /scratch

## Step 2: clone repository from github
    git clone https://github.com/redtreevole/Foss2020Team6project

## Step 3: stage data in from Data Store at iplant/home/filiaultd/GWAS_data to directory 001.data
    mkdir 001.data
    cd 001.data
    irsync -r i:GWAS_data ./
 
 Data files required are:
 * genotypes in bed format (with accompanying bim and fam files) - __02_2.3M_200Swedes.biallelic.imputed.filtered.bed__
 * phenotype TXT file with accessions in same order as genotypes, phenotype values given in columns with appropriate column headers - __{pheno_name}.fitness.txt__
 * K matrix calculated from genotype data - __K.matrix.200Swedes.labels.txt__

## Step 3: run snakemake pipeline
    cd ..
    snakemake --use-conda --cores 1
### I can't make the conda environment!!! :confounded:

## Step 4: stageout results to Data Store
    cd 003.results
    irsync -r /scratch/Foss2020Team6project/limixGWAS/003.results i:GWAS_data/003.results

## Other general notes
use cyberduck to move data to/from local computer to Data Store
