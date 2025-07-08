#!/bin/bash

## LVRS GWAS on dentition phenotypes (four)
## pooja.singh09@gmail.com
## dec 2024 


# 2_gwasdentition.sbatch
	#SBATCH -J "2_gwascol"
	#SBATCH --mem-per-cpu=8000
	#SBATCH --time=24:00:00
	#SBATCH -n 1
	#SBATCH -o 2_gwascol.log
	#SBATCH -e 2_gwascol.err

	# Change working directory to global scratch
	cd /cluster/scratch/posingh/gwas/gemma_dentition
	source /cluster/project/gdc/shared/stack/GDCstack.sh
	module load module load gemma/0.98.5

	#snppath=/cluster/work/gdc/shared/p929/victoriaGenomes/vcf_pnye3/2023_old_without_GQ/
	gwaspath=/cluster/scratch/posingh/gwas/gemma_dentition

	# Replace plink .fam file with phenotypes (GEMMA only reads second and 6th column, and 6th column onwards should be the phenotype
	
	#sed 1d $gwaspath/gwas_dentition_phenotypes_v1.lvrs.NA  | awk '{print $1,$1,0,0,0,$2,$3,$4,$5}' > $gwaspath/all.SNPs.filt2.dentition.imp.fam

	# Run GEMMA: first compute relatedness matrix

	gemma -bfile $gwaspath/all.SNPs.filt2.dentition.imp -gk 1 -o all.SNPs.filt2.dentition.imp
posingh@eu-login-12:/cluster/scratch/posingh/gwas/gemma_dentition$ cat 3_gwasdentition_10runs.sh

max=10
for i in `seq 1 $max`
do
    echo "sbatch -n 1 --mem-per-cpu=30000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt -n 1 -bslmm 1 -w 5000000 -s 20000000 -o all.SNPs.filt.dentition.imp.bslmm.dentition.CuspShape.$i'"

echo "sbatch -n 1 --mem-per-cpu=30000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt -n 2 -bslmm 1 -w 5000000 -s 20000000 -o all.SNPs.filt.dentition.imp.bslmm.dentition.InnertoothrowUJ.$i'"

echo "sbatch -n 1 --mem-per-cpu=30000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt -n 3 -bslmm 1 -w 5000000 -s 20000000 -o all.SNPs.filt.dentition.imp.bslmm.dentition.LPJtoothMolarisation.$i'"

echo "sbatch -n 1 --mem-per-cpu=30000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt -n 4 -bslmm 1 -w 5000000 -s 20000000 -o all.SNPs.filt.dentition.imp.bslmm.dentition.LPJkeelDepth.$i'"


done

####
# run output file using sbatch -n 1 --mem-per-cpu=30000 --time=12:00:00 -
#
