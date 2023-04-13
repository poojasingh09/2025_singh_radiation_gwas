## LVRS GWAS on body colour and patterning and head shape PCs from Kass
## pooja.singh09@gmail.com

# 1_gwas_jaw.lsf
	#BSUB -J "1_gwasjaw"
	#BSUB -R "rusage[mem=2500]"
	#BSUB -W 24:00
	#BSUB -n 4
	#BSUB -o 1gwasjaw.log
	#BSUB -e 1gwasjaw.err

	# Change working directory to global scratch
	cd /cluster/scratch/posingh/gwas/jawshape

	module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.8 java/1.8.0_73 beagle/4.1 gemma/0.94 vcftools/0.1.16 zlib/1.2.8 openblas/0.2.13_seq plink/1.90
	perl-init

  #set paths
  
	snppath=/cluster/scratch/posingh/bcftools/v3
	gwaspath=/cluster/scratch/posingh/gwas/jawshape

  # get list of indivduals
  
	fgrep -w -f lvrs_samples.txt gwas_headshape_PCall_edited.txt > gwas_headshape_PCall_edited.lvrs
	cut -f 1 gwas_headshape_PCall_edited.lvrs > gwas_headshape_PCall_edited.lvrs.inds

	# 1) Extract subset of individuals and filter VCF file for SNPs with <10% missing dat

	bcftools view --force-samples -S $gwaspath/gwas_headshape_PCall_edited.lvrs.inds \
	$snppath/all.SNPs.raw.vcf.gz | bgzip -c > $gwaspath/all.SNPs.raw.inds.vcf.gz 
	

	#10 % missing indiviuals only, biallic
	bcftools +fill-tags | bcftools view -i 'MAF>0.05 & AN/(2*N_SAMPLES)>0.9 & N_ALT=1 & FORMAT/DP > 20 & INFO/DP > 20 & QUAL > 30' $gwaspath/all.SNPs.raw.inds.vcf.gz | bgzip -c > $gwaspath/all.SNPs.filt.jaw.vcf.gz
	tabix -p vcf $gwaspath/all.SNPs.filt.jaw.vcf.gz

	# 2) Impute missing genotypes with BEAGLE
	beagle gt=$gwaspath/all.SNPs.filt.jaw.vcf.gz impute=T nthreads=4 out=$gwaspath/all.SNPs.filt.jaw.imp
	
	tabix -p vcf $gwaspath/all.SNPs.filt.jaw.imp.vcf.gz
	

	# 3 	# Combine files and convert genotypes to binary PED format # First from VCF to PLINK format
	plink --vcf $gwaspath/all.SNPs.filt.jaw.imp.vcf.gz --make-bed --out $gwaspath/all.SNPs.filt.jaw.imp
