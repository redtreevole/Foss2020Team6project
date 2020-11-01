phenos=["ADA_2011","ADA_2012"]
#, "RAT_2011", "ULL_2011", "RAM_2012",  "ADA_2012", "RAT_2012", "ULL_2012"]

rule all:
	input:
		expand("003.results/{pheno}.limix.results.csv", pheno=phenos)

rule limix:
	input:
		pheno_file="001.data/{pheno}.fitness.txt",
		geno_file="001.data/02_2.3M_200Swedes.biallelic.imputed.filtered.bed",
		kin_file="001.data/K.matrix.200Swedes.labels.txt",
		#pheno_name=PHENOS
	output:
		#expand("/003.results/{pheno_name}.limix.results.csv", pheno_name=PHENOS)
		"003.results/{pheno}.limix.results.csv"		
	conda:
        	"limix.environ.yml"	
	shell:
		"python limix.marginal.py -pf {input.pheno_file} -gf {input.geno_file} -kf {input.kin_file} -pn {wildcards.pheno}"
