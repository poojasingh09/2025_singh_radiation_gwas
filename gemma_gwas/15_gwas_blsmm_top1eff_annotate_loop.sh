for j in `cat trait`;

do

Rscript --vanilla 14_gwas_blsmm_top1eff_annotate.R $j

done
