# bsub < 2_gwasjaw.lsf
#segmentation fault error may arise if individuals in vcf dont match individuals in phenotype list

	#BSUB -J "gwas2jaw"
	#BSUB -R "rusage[mem=4000]"
	#BSUB -W 4:00
	#BSUB -n 1
	#BSUB -o 2_gwasjaw.log
	#BSUB -e 2_gwasjaw.err

	module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.8 java/1.8.0_73 beagle/4.1 gemma/0.94 vcftools/0.1.16 zlib/1.2.8 openblas/0.2.13_seq plink/1.90
	perl-init

	gwaspath=/cluster/scratch/posingh/gwas/jawshape

	# Replace plink .fam file with phenotypes (GEMMA only reads second and 6th column, and 6th column onwards should be the phenotype
	#sed 1d $gwaspath/gwas_headshape_PCall_edited.txt  | awk '{print $1,$1,0,0,0,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' > $gwaspath/all.SNPs.filt.jaw.imp.fam
	
	#fgrep -w -f $gwaspath/gwas_headshape_PCall_edited.lvrs.inds $gwaspath/all.SNPs.filt.jaw.imp.fam > $gwaspath/a.txt
	#mv $gwaspath/a.txt $gwaspath/all.SNPs.filt.jaw.imp.fam

	# Run GEMMA: first compute relatedness matrix

	gemma -bfile $gwaspath/all.SNPs.filt.jaw.imp -gk 1 -o all.SNPs.filt.jaw.imp
  
  
  # RUN GEMMA:
  
  # bsub < 3_gwas_jaw.lsf

	#BSUB -J "gwas3a[1-2]%2"
	#BSUB -R "rusage[mem=4000]"
	#BSUB -W 180:00
	#BSUB -n 1
	#BSUB -o 3_gwasjaw.log.%J.%I
	#BSUB -e 3_gwasjaw.err.%J.%I
	
	module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.8 java/1.8.0_73 beagle/4.1 gemma/0.94 vcftools/0.1.16 zlib/1.2.8 openblas/0.2.13_seq plink/1.90
	perl-init

	cd /cluster/scratch/posingh/gwas/jawshape

	# Run GEMMA: lmm for all traits (ie head shape PCs) separately

	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 1 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC1
  gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 2 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC2
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 3 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC3
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 4 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC4
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 5 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC5
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 6 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC6
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 7 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC7
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 8 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC8
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 9 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC9
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 10 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC10
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 11 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC11
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 12 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC12
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 13 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC13
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 14 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC14
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 15 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC15
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 16 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC16
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 17 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC17
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 18 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC18
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 19 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC19
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt  -n 20 -lmm 4 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC20


        

	# Run GEMAA: mvlmm multivariate lmm for multiple traits (ie PCs) at once
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt -lmm 4 -n 1 2 -o all.SNPs.filt.jaw.imp.lmm.jaw.PC1_2
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt -lmm 4 -n 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 -o all.SNPs.filt.jaw.imp.lmm.jaw.PCall
