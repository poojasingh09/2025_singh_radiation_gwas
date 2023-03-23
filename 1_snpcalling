# pooja.singh09@gmail.com
# 2023
# bcftools_callsnps.lsf

#BSUB -J "callg[1-22]%22"
#BSUB -R "rusage[mem=6000]"
#BSUB -W 120:00
#BSUB -n 1
#BSUB -o snpcallg.log.%I.%J
#BSUB -e snpcallg.err.%I.%J

cd /cluster/scratch/posingh/bcftools/v3
path=/cluster/scratch/posingh/bcftools/v3

module load gcc/4.8.2 gdc samtools/1.8 vcftools/0.1.16

chr=$LSB_JOBINDEX

### run bcftool on all bam files cohort using P.nye v3 reference genome

bcftools mpileup -Ou -I -r chr${chr} -f /cluster/work/gdc/shared/p929/ref-genome/PunNye3.0.fasta /cluster/work/gdc/shared/p929/victoriaGenomes/bamFilesPnye3/*bam \
-a "DP,AD" | bcftools call -Ob -mv -f GQ > ${path}/all.${chr}.raw.bcf.gz
