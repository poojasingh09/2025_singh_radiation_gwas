 bsub < 4_gwas_jaw.lsf

	#BSUB -J "gwas4a[1-200]%20"
	#BSUB -R "rusage[mem=8000]"
	#BSUB -W 12:00
	#BSUB -n 1
	#BSUB -o 4gwasjaw.log.%J.%I
	#BSUB -e 4gwasjaw.err.%J.%I
	
	module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.8 java/1.8.0_73 beagle/4.1 gemma/0.94 vcftools/0.1.16 zlib/1.2.8 openblas/0.2.13_seq plink/1.90
	perl-init

	cd /cluster/scratch/posingh/gwas/jawshape

	#head -n1 gwas_headshape_PC12_edited.lvrs | sed 's/\t/\n/g' | sed 1d | perl -ne 'for$i(0..9){print}' | awk  'FNR==0{print $0;next}p!=$1{c=0;p=$1}{print $0"."++c}' > jawloop
	
	IDX=$LSB_JOBINDEX
	trait=`sed -n ${IDX}p <jawloop`
	col=`sed -n ${IDX}p <jawloopcol` #jawloopcol is a file with one column with  numbers 1 and 2 (because there are 2 traits, PC1 and PC2) each number repeated 10 times. 


	# Run Gemma: bayesian sparse linear mixed model (bslmm) for modelling polygenic traits with simple and complex architectures
	gemma -bfile all.SNPs.filt.jaw.imp -k ./output/all.SNPs.filt.jaw.imp.cXX.txt -n ${col} -bslmm 1 -o all.SNPs.filt.jaw.imp.bslmm.jaw.${trait}
