#BSUB -J "ant-fst"
#BSUB -R "rusage[mem=1000]"
#BSUB -W 1:00
#BSUB -n 1


cd /cluster/scratch/posingh/
module load gcc/4.8.2 gdc samtools/1.8 vcftools/0.1.16

vcftools --vcf /cichlid.vcf --fst-window-size 20000 --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --out pop1_pop2.fst.out 2>&1 | tee -a   pop1_pop2.fst.log

