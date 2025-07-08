#!/bin/bash
## LVRS GWAS phylogeny
## pooja.singh09@gmail.com
## April 2025


#SBATCH -J "2_filter_vcf_parallel"
#SBATCH --array=1-22%22
#SBATCH --mem-per-cpu=4000
#SBATCH --time=24:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --tmp=120G
#SBATCH --output=2_filter_vcf_parallel_%a.log
#SBATCH --err=2_filter_vcf_parallel_%a.err

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load stack/.2024-06-silent  gcc/12.2.0 bcftools/1.20 vcftools/0.1.16-kpfcvxs samtools/1.20

cd /cluster/scratch/posingh/gwas_phylogeny


IDX=$SLURM_ARRAY_TASK_ID

#bcftools view all.raw.bcf.gz -Oz -o all.raw.vcf.gz
vcftools --gzvcf all.raw.vcf.gz --chr chr${IDX} --keep all_bam_wanted_samples_fishec --minGQ 30 --minDP 20 --min-meanDP 20 --maf 0.05 --max-missing 0.8 --max-alleles 2 --remove-indels --recode --stdout | bgzip -c > all.SNPs.${IDX}.filt.vcf.gz


# concat all chrs
bcftools concat all.SNPs.1.filt.vcf.gz all.SNPs.2.filt.vcf.gz all.SNPs.3.filt.vcf.gz all.SNPs.4.filt.vcf.gz all.SNPs.5.filt.vcf.gz all.SNPs.6.filt.vcf.gz all.SNPs.7.filt.vcf.gz all.SNPs.8.filt.vcf.gz all.SNPs.9.filt.vcf.gz all.SNPs.10.filt.vcf.gz all.SNPs.11.filt.vcf.gz all.SNPs.12.filt.vcf.gz all.SNPs.13.filt.vcf.gz all.SNPs.14.filt.vcf.gz all.SNPs.15.filt.vcf.gz all.SNPs.16.filt.vcf.gz all.SNPs.17.filt.vcf.gz all.SNPs.18.filt.vcf.gz all.SNPs.19.filt.vcf.gz all.SNPs.20.filt.vcf.gz all.SNPs.21.filt.vcf.gz all.SNPs.22.filt.vcf.gz -Oz -o all.filt.vcf.gz
##################################################################################################################################################################################################
######control header
echo "Starting ${SLURM_JOB_ID} at $(date)"
echo "Pooja Singh"
##################################################################################################################################################################################################
#### 09.04.25



file=all.filt.vcf.gz

echo "Start processing ${file}"

## keep sites where at least one individual is 00 and 11, nonvariant sites are uninformative

bcftools view -i 'COUNT(GT="RR")>0 & COUNT(GT="AA")>0'  $file  -Oz -o all.filt.vcf.var.gz
tabix -fp vcf all.filt.vcf.var.gz

## thin SNPs, this speeds up phylogenetic reconstruction by removing linked SNPs that are unlikely to provide extra information for phylogenetic reconstruction

vcftools --gzvcf all.filt.vcf.var.gz \
 --thin 1000 --stdout --recode | gzip > all.filt.vcf.var.1kbthinned.gz


## vcf to alignment in phylip format
#git clone https://github.com/edgardomortiz/vcf2phylip.git # this command you can run without bsub as it just downloads the vcf2phylip script from github

./vcf2phylip/vcf2phylip.py -i all.filt.vcf.var.1kbthinned.gz -o all.filt.vcf.var.1kbthinned.phy

##iqtreeee

iqtree2 -s all.filt.vcf.var.1kbthinned.phy -m GTR+G -bb 1000 -alrt 1000 -pre all.filt.vcf.var.1kbthinned

##################################################################################################################################################################################################
####job controlling, grep "Job:" *.out
seff $SLURM_JOB_ID

echo "Job: ${SLURM_JOB_ID} successfully finished $(date)"
# happy end
exit 0
##################################################################################################################################################################################################



