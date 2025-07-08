# pooja.singh09@gmail.com
# dec 2024
# in case you also want to run lmm

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load gemma/0.98.5

	# Run GEMMA: lmm

	sed 1d gwas_dentition_phenotypes_v1.lvrs.NA  | awk '{print $1,$1,0,0,0,$2,$3,$4,$5}' > all.SNPs.filt2.dentition.imp.fam
	sbatch -n 1 --mem-per-cpu=10000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt  -n 1 -lmm 4 -o all.SNPs.filt2.dentition.imp.lmm.dentition.CuspShape'
        
	sed 1d gwas_dentition_phenotypes_v1.lvrs.NA  | awk '{print $1,$1,0,0,0,$3,$2,$4,$5}' > all.SNPs.filt2.dentition.imp.fam
	sbatch -n 1 --mem-per-cpu=10000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt  -n 1 -lmm 4 -o all.SNPs.filt2.dentition.imp.lmm.dentition.InnertoothrowUJ'

	sed 1d gwas_dentition_phenotypes_v1.lvrs.NA  | awk '{print $1,$1,0,0,0,$4,$2,$3,$5}' > all.SNPs.filt2.dentition.imp.fam
	sbatch -n 1 --mem-per-cpu=10000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt  -n 1 -lmm 4 -o all.SNPs.filt2.dentition.imp.lmm.dentition.LPJtoothMolarisation'

	sed 1d gwas_dentition_phenotypes_v1.lvrs.NA  | awk '{print $1,$1,0,0,0,$5,$2,$3,$4}' > all.SNPs.filt2.dentition.imp.fam
	sbatch -n 1 --mem-per-cpu=10000 --time=12:00:00 --wrap='gemma -bfile all.SNPs.filt2.dentition.imp -k ./output/all.SNPs.filt2.dentition.imp.cXX.txt  -n 1 -lmm 4 -o all.SNPs.filt2.dentition.imp.lmm.dentition.LPJkeelDepth'
