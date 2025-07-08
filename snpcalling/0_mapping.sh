# usage bsub < bwaalign.lsf
# bwaalign.lsf
	#!/bin/bash
	#BSUB -P mapping
	#BSUB -J "bwaalign[1-3]%3"
	#BSUB -R "rusage[mem=1000]"
	#BSUB -n 4
	#BSUB -W 48:00
	#BSUB -o log.bwaalign.%J.%I
	#BSUB -e err.bwaalign.%J.%I
	
	# Get individual number from file

	#cp /cluster/work/gdc/shared/p500/victoriaGenomes/fastq/2022/list.txt .

	IDX=$LSB_JOBINDEX

	ind=`sed -n ${IDX}p <list.txt`

	# Load the required modules
	module load gcc/4.8.2 gdc bwa/0.7.17 perl/5.18.4 samtools/1.9 fastp/0.20.0 fastqc/0.11.4
	perl-init

	# Define paths to reference genome, fastq file folder and alignment destination folder
	ref="/cluster/work/gdc/shared/p929/ref-genome/pnye3.0"
	#fastqpath="/cluster/work/gdc/shared/p500/victoriaGenomes/fastq/2022/"
	#bampath="/cluster/work/gdc/shared/p500/victoriaGenomes/bamFilesPnye3/2022/"
	outpath="/cluster/scratch/posingh/reseq_comparison/"	

	# Change working directory to local scratch
	cd /scratch

	# Run BWA MEM alignment and stream output to bam format
	#    Options:
	#    -t 8, use of 8 threads, hyper-threading on Euler activated, requesting 4 cores, but using 8 threads
	#    -T 20, don't output alignments with MAPQ < 20
	#    -R '', append read group information to header
	#  samtools: -F 256, filter out secondary alignments
	bwa mem -t 4 -T 20 \
	-R '@RG\tSM:'${ind}'\tID:'${ind}'\tPL:ILLUMINA\tPG:bwa-mem_0.7.17' \
	${ref} ${outpath}${ind}_R1.trimmed.fastq.gz ${outpath}${ind}_R2.trimmed.fastq.gz 2> ${outpath}${ind}.bwamem.log | \
	samtools view -bF 256 > ${outpath}${ind}.bwamem.bam

	# Sort and index bam file
	#samtools sort -O BAM -T ./${ind}.bwamem ${ind}.bwamem.bam > ${ind}.bwamem.sorted.bam
	#mv ${ind}.bwamem.sorted.bam ${ind}.bwamem.bam
	#samtools index ${ind}.bwamem.bam

	# Upload bam file, index and logfile to alignment destination folder
	#mv ${ind}.bwamem.* ${bampath}
