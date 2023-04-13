# pooja.singh09@gmail.com
# filt for baseline parameters. PLEASE FILTER MORE FOR YOUR SPECIFIC DOWNSTREAM ANALYSES!!!!! ie gwas, demographics etc

#BSUB -J "snpfilt[1-22]%22" 
#BSUB -W 24:00 
#BSUB -R "rusage[mem=4000]" 
#BSUB -o snpfilt1.log.%J.%I
#BSUB -e snpfilt1.err.%J.%I


#load modules
module load gcc/4.8.2 gdc samtools/1.8 vcftools/0.1.16

#set chr names
chr=$LSB_JOBINDEX

#set file name prefix
prefix=all.${chr}.raw

#change dir
cd /cluster/scratch/posingh/bcftools/v3

## SNPs only
vcftools --bcf $prefix.bcf.gz --remove-indels --recode --recode-INFO-all --stdout |  bgzip -c >  $prefix.snps.vcf.gz
tabix $prefix.snps.vcf.gz

## depth, qual, maf filters
vcftools --gzvcf all.SNPs.raw.vcf.gz --minGQ 30 --minDP 20 --min-meanDP 20 --maf 0.05 --max-missing 0.8 --max-alleles 2 --recode-INFO-all --recode --stdout | bgzip -c > all.SNPs.filt.vcf.gz
