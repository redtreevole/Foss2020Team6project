#!/usr/bin/env python
# coding: utf-8

# # Introduction
# Script to run marginal GWAS in limix with field experiment SNPs 
# DLF 11 June 20
# reworked from ipynb to .py script 01Sept20

#########################
### 0. Set up environment

from limix.qtl import scan
import pandas as pd
import numpy as np
import h5py
from bisect import bisect
from limix import plot
import random
from sklearn.decomposition import PCA
import argparse
import os
from statsmodels.stats import multitest
import math
import limix


##############################
### 1. parse arguments
parser = argparse.ArgumentParser(description = 'Parse parameters for marginal GWAS')
parser.add_argument('-pf', '--phenfile', help = 'Path to phenotype file. TXT file with accessions in same order as genotypes, phenotype values given in columns with appropriate column headers. Filename is used as phenotype name.', required = True)
parser.add_argument('-gf', '--genfile', help = 'Path to genotype file in bed format', required = True)
parser.add_argument('-kf', '--kfile', help = 'Path to K matrix file', required=True)
parser.add_argument('-pn', '--phenname', help = 'Specify the column name of the phenotype in the phenotype file', required = True)
args = parser.parse_args()


##############################
### 2.  Settings

### NOTE: this script assumes that accessions are in the same order in all these files!
### There is a check of this in the script, but it might still run if this condition is not met.
### So check this beforehand (or at least check the log file to make sure condition is met)

# phenotype
pheno_file = args.phenfile
#pheno_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/slopeNS.blups.txt'

# genotype
geno_file = args.genfile
#geno_file = '/groups/nordborg/projects/field_experiments/001.common.reference.files/006.field.SNPs.Fernando/imputed/02_2.3M_200Swedes.biallelic.imputed.filtered.bed'

# kinship matrix
kin_file = args.kfile
#kin_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/K.matrix.200Swedes.labels.txt'
# minor allele frequency threshold, not set in this script.  This was set when generating the input bim file
#MAF_thrs = 0.1 

# pheno name and output results
#pheno_name = "mu2.sl"
pheno_name = args.phenname
output_file = '/groups/nordborg/projects/field_experiments/adaptation_sweden/common.gardens/data/limix/' + pheno_name +'.limix.results.csv'
# setting limix file as default output location.  Change to add this as a proper argument if needed/desired


#############################
### 3. Read phenotypes

# Read in phenotype file and subset by phenotype of interest
phenoAll = pd.read_csv(pheno_file, index_col = 1, sep=" ")
phenoAll.index = phenoAll.index.map(str) ## need to make index type match that of genotype index so can check if in same order

pheno = phenoAll[['id',pheno_name]]
# get indices of NA values to filter other data later
ind = np.where(pheno[pheno_name].isnull())

# remove NA values
pheno = pheno[np.isfinite(pheno)]
pheno = pheno.dropna()
pheno_ids = pheno[['id']]
pheno = pheno[[pheno_name]]
#print(pheno.shape)
#print(pheno.head(n=5))

#print(pheno.isnull().any().any()) #checks for NAs
Y=pheno.to_numpy()
print(Y.shape)


####################################
### 4. Read genotypes

import limix
(bim, fam, bed) = limix.io.plink.read(geno_file, verbose=False)
G=bed.compute()
G=G.transpose()
#remove NA accessions noted with ind giving rownumber
G=np.delete(G, ind, axis=0)


####################################
### 5. Read K matrix

K = pd.read_csv(kin_file, index_col = 0, sep=" ")
K.index = K.index.map(str)
Kp = K
K = K.to_numpy()

#remove any NA accessions with index
#note, this only works because I already know that all these matrices are in the same order
#if writing for different data or troubleshooting, might need to address this
K=np.delete(K, ind, axis=0)
K=np.delete(K, ind, axis=1)


#########################################
### 6. Make sure all datasets in same accession order

### THIS SHOULD BE CHECKED BEFOREHAND!!!!

# also remove NA row from fam, which gives information about genotypes
fam.drop(fam.index[ind], inplace=True)
# and from Kp, which gives K accession info
Kp.drop(Kp.index[ind], inplace=True)
Kp.drop(Kp.columns[ind], axis=1, inplace=True)


GP_check = fam.index == pheno_ids.index
PK_check = pheno_ids.index == Kp.index
if np.count_nonzero(GP_check == False) != 0 and np.count_nonzero(PK_check == False) != 0 : print("CAUTION: Not all data files are in the same order!")
if np.count_nonzero(GP_check == False) == 0 and np.count_nonzero(PK_check == False) == 0 : print("All input files in same order")

    
######################################
### 7. Run marginal GWAS in limix

r = scan(G=G, Y=Y, K = K, lik = 'normal', M = None, verbose = True)

####################################
### 8. Output results

chrom = bim[['chrom']]
pos = bim[['pos']]
#extract pvals
pvalues = r.stats.pv20.tolist()
#extract effect sizes
effsizes = r.effsizes['h2']['effsize'][r.effsizes['h2']['effect_type'] == 'candidate'].tolist()

gwas_results = np.c_[chrom, pos, pvalues, effsizes]
gwas_results = pd.DataFrame(data=gwas_results, index=None, columns=["chrom", "pos", "pv", "GVE"])
gwas_results.to_csv(output_file, index = False)

########################################
### 9. Manhattan plot

#gwas_results.dtypes
gwas_results["chrom"] = gwas_results['chrom'].astype('int')
gwas_results["pos"] = gwas_results['pos'].astype('int')
gwas_results["pv"] = gwas_results['pv'].astype('float')
gwas_results["GVE"] = gwas_results['GVE'].astype('float')
#gwas_results.dtypes

Bonferroni = multitest.multipletests(pvalues, alpha = 0.05, method = 'fdr_bh')[3]
plot.manhattan(gwas_results)
plt = plot.get_pyplot() 
_ = plt.axhline(-math.log(Bonferroni, 10), color='red')  
plt.savefig(f'{output_file}.manhattan.png')
plt.close()

###########################################
### 10. QQplot

plot.qqplot(gwas_results.pv)
plt = plot.get_pyplot()
plt.savefig(f'{output_file}.QQplot.png')
plt.close()


###########################################
### 11. Calculate heritability

from numpy.random import RandomState
from limix.her import estimate

herit = estimate(y=Y, K=K, lik="normal", verbose=False)
print(herit)
f = open(f"{output_file}.herit.txt", "w")
f.write(str(herit))
f.close()

